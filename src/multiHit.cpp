#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <random>
#include <set>
#include <vector>

#include "multiHit.h"
#include "utils.h"

#ifdef ENABLE_PROFILE
static double iter_base_times[TIMING_COUNT];

static inline void snapshot_iter_timers() {
  for (int i = 0; i < TIMING_COUNT; ++i)
    iter_base_times[i] = elapsed_times[i];
}

static inline double iter_delta(int timer_id) {
  double d = elapsed_times[timer_id] - iter_base_times[timer_id];
  return d < 0.0 ? 0.0 : d;
}
#endif

//////////////////////////////  Start Allreduce_hierarchical
/////////////////////////

#ifdef HIERARCHICAL_COMMS
#include <unistd.h>

#define ALL_REDUCE_FUNC Allreduce_hierarchical
#define EXECUTE execute_hierarchical
static void Allreduce_hierarchical(void *sendbuf, void *recvbuf, int count,
                                   MPI_Datatype datatype, MPI_Op op,
                                   CommsStruct &comms) {

  // Datatype size
  int datatype_size;
  MPI_Type_size(datatype, &datatype_size);

  // Allocate buffer for local reduction
  void *local_result = malloc(count * datatype_size);

  START_TIMING(local_reduce);
  // Phase 1: Local Reduction
  MPI_Reduce(sendbuf, local_result, count, datatype, op, 0, comms.local_comm);
  END_TIMING(local_reduce, elapsed_times[COMM_LOCAL_TIME]);
  // Allocate buffer for global reduction (only needed for rank 0 in
  // local_comm)
  void *global_result = NULL;
  if (comms.is_leader) {
    global_result = malloc(count * datatype_size);
  }

  // Phase 2: Global Reduction (only rank 0 in each node participates)
  if (comms.is_leader) {
    START_TIMING(global_reduce);
    MPI_Allreduce(local_result, global_result, count, datatype, op,
                  comms.global_comm);
    END_TIMING(global_reduce, elapsed_times[COMM_GLOBAL_TIME]);
    // Copy the final result to recvbuf
    memcpy(recvbuf, global_result, count * datatype_size);
  }

  // Phase 3: Broadcast result to all processes in the local communicator
  START_TIMING(local_bcast);
  MPI_Bcast(recvbuf, count, datatype, 0, comms.local_comm);
  END_TIMING(local_bcast, elapsed_times[COMM_LOCAL_TIME]);

  // Cleanup
  free(local_result);
  if (comms.is_leader)
    free(global_result);
}

static inline WorkChunk calculate_node_range(LAMBDA_TYPE num_Comb,
                                             const CommsStruct &comms) {
  const LAMBDA_TYPE base = num_Comb / comms.num_nodes;
  const LAMBDA_TYPE extra = num_Comb % comms.num_nodes;
  const int k = comms.my_node_id;

  const LAMBDA_TYPE start = k * base + std::min<LAMBDA_TYPE>(k, extra);
  const LAMBDA_TYPE len = base + (k < extra ? 1 : 0);

  return {start, start + len - 1};
}

static inline WorkChunk calculate_worker_range(const WorkChunk &leaderRange,
                                               int worker_id, int num_workers) {

  const LAMBDA_TYPE start = leaderRange.start + worker_id * CHUNK_SIZE;

  if (start > leaderRange.end) {
    return {0, -1};
  }

  const LAMBDA_TYPE end = std::min(start + CHUNK_SIZE - 1, leaderRange.end);

  return {start, end};
}

inline LAMBDA_TYPE length(const WorkChunk &c) { return c.end - c.start + 1; }

inline static void handle_local_work_steal(WorkChunk &availableWork,
                                           MPI_Status st,
                                           LAMBDA_TYPE &combs_dispensed,
                                           const CommsStruct &comms) {
  int requester = st.MPI_SOURCE;
  char dummy;
  START_TIMING(recv);
  MPI_Recv(&dummy, 1, MPI_BYTE, requester, TAG_REQUEST_WORK, comms.local_comm,
           MPI_STATUS_IGNORE);
  END_TIMING(recv, elapsed_times[COMM_LOCAL_TIME]);
  WorkChunk reply;
  if (availableWork.start <= availableWork.end) {
    LAMBDA_TYPE chunkEnd =
        std::min(availableWork.start + CHUNK_SIZE - 1, availableWork.end);
    reply = {availableWork.start, chunkEnd};
    availableWork.start = (chunkEnd + 1);
    combs_dispensed += length(reply);

  } else {
    reply = {0, -1};
  }

  START_TIMING(send);
  MPI_Send(&reply, sizeof(WorkChunk), MPI_BYTE, requester, TAG_ASSIGN_WORK,
           comms.local_comm);
  END_TIMING(send, elapsed_times[COMM_LOCAL_TIME]);
}

inline static void inter_node_work_steal_victim(WorkChunk &availableWork,
                                                MPI_Status st, int &my_color,
                                                Token &tok,
                                                const CommsStruct &comms) {
  char dummy;
  MPI_Request rq_recv;
  START_TIMING(victim_irecv);
  MPI_Irecv(&dummy, 1, MPI_BYTE, st.MPI_SOURCE, TAG_NODE_STEAL_REQ,
            comms.global_comm, &rq_recv);
  END_TIMING(victim_irecv, elapsed_times[COMM_GLOBAL_TIME]);

  WorkChunk reply;

  if (availableWork.start < availableWork.end) {
    LAMBDA_TYPE len = availableWork.end - availableWork.start + 1;
    LAMBDA_TYPE mid = availableWork.start + len / 2;
    reply = {mid, availableWork.end};
    availableWork.end = mid - 1;
    my_color = BLACK;
  } else {
    reply = {0, -1};
  }

  MPI_Request rq;
  START_TIMING(victim_isend);
  MPI_Isend(&reply, sizeof(WorkChunk), MPI_BYTE, st.MPI_SOURCE,
            TAG_NODE_STEAL_REPLY, comms.global_comm, &rq);
  END_TIMING(victim_isend, elapsed_times[COMM_GLOBAL_TIME]);
}

static inline void root_broadcast_termination(const CommsStruct &comms,
                                              MPI_Win &term_win) {
  bool termination_signal = true;
  for (int rank = 0; rank < comms.num_nodes; ++rank) {
    START_TIMING(node_term);
    MPI_Put(&termination_signal, 1, MPI_C_BOOL, rank, 0, 1, MPI_C_BOOL,
            term_win);
    END_TIMING(node_term, elapsed_times[COMM_GLOBAL_TIME]);
  }
  MPI_Win_flush_all(term_win);
}

static inline void try_forward_token_if_idle(WorkChunk &availableWork,
                                             bool &have_token, int &my_color,
                                             Token &tok, const int next_leader,
                                             MPI_Win &term_win,
                                             const CommsStruct &comms) {

  if (!have_token || length(availableWork) > 0) {
    return;
  }

  if (my_color == BLACK)
    tok.colour = BLACK;

  if (comms.global_rank == 0) {
    if (tok.colour == WHITE && tok.finalRound) {
      root_broadcast_termination(comms, term_win);
    } else {
      tok.finalRound = (tok.colour == WHITE);
      tok.colour = WHITE;
      START_TIMING(send);
      MPI_Send(&tok, sizeof(Token), MPI_BYTE, next_leader, TAG_TOKEN,
               comms.global_comm);
      END_TIMING(send, elapsed_times[COMM_GLOBAL_TIME]);
    }
  } else {
    START_TIMING(send);
    MPI_Send(&tok, sizeof(Token), MPI_BYTE, next_leader, TAG_TOKEN,
             comms.global_comm);
    END_TIMING(send, elapsed_times[COMM_GLOBAL_TIME]);
  }

  have_token = false;
  my_color = WHITE;
}

inline static void receive_token(Token &tok, MPI_Status st, bool &have_token,
                                 const CommsStruct &comms) {
  MPI_Status status;
  START_TIMING(recv);
  MPI_Recv(&tok, sizeof(Token), MPI_BYTE, st.MPI_SOURCE, TAG_TOKEN,
           comms.global_comm, &status);
  END_TIMING(recv, elapsed_times[COMM_GLOBAL_TIME]);
  have_token = true;
}

inline static WorkChunk
assign_and_update_availableWork(const WorkChunk &loot, int num_workers,
                                LAMBDA_TYPE &combs_dispensed,
                                const CommsStruct &comms) {
  LAMBDA_TYPE max_end = loot.start - 1;
  for (int w = 1; w <= num_workers; ++w) {
    WorkChunk work = calculate_worker_range(loot, w - 1, num_workers);
    START_TIMING(assign_work);
    MPI_Send(&work, sizeof(WorkChunk), MPI_BYTE, w, TAG_ASSIGN_WORK,
             comms.local_comm);
    END_TIMING(assign_work, elapsed_times[COMM_LOCAL_TIME]);
    combs_dispensed += length(work);
    if (work.end > max_end)
      max_end = work.end;
  }
  if (max_end >= loot.end) {
    return {0, -1};
  } else {
    return {max_end + 1, loot.end};
  }
}

inline static void inter_node_work_steal_initiate(
    WorkChunk &availableWork, MPI_Status st, int num_workers, int &my_color,
    Token &tok, MPI_Win &term_win, bool *global_done,
    LAMBDA_TYPE &combs_dispensed, const CommsStruct &comms) {
  int myRank = comms.global_rank;
  int nLeaders = comms.num_nodes;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> gen(0, nLeaders - 1);
  bool lootReceived = false;
  char dummy;

  for (int attempt = 0; attempt < NUM_RETRIES && !lootReceived; ++attempt) {
    int victim;
    do {
      victim = gen(rng);
    } while (victim == myRank);

    MPI_Request rq;
    MPI_Status status;
    START_TIMING(isend);
    MPI_Isend(&dummy, 1, MPI_BYTE, victim, TAG_NODE_STEAL_REQ,
              comms.global_comm, &rq);
    END_TIMING(isend, elapsed_times[COMM_GLOBAL_TIME]);

    WorkChunk loot;
    MPI_Request rq_recv;
    START_TIMING(irecv);
    MPI_Irecv(&loot, sizeof(WorkChunk), MPI_BYTE, victim, TAG_NODE_STEAL_REPLY,
              comms.global_comm, &rq_recv);
    END_TIMING(irecv, elapsed_times[COMM_GLOBAL_TIME]);

    int completed = 0;
    MPI_Win_sync(term_win);

    while (!completed && !(*global_done)) {
      START_TIMING(recv_test);
      MPI_Test(&rq_recv, &completed, MPI_STATUS_IGNORE);
      END_TIMING(recv_test, elapsed_times[COMM_GLOBAL_TIME]);
      int flag = 0;
      MPI_Status st;
      START_TIMING(iprobe);
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comms.global_comm, &flag, &st);
      END_TIMING(iprobe, elapsed_times[COMM_GLOBAL_TIME]);

      if (flag) {
        switch (st.MPI_TAG) {

        case TAG_NODE_STEAL_REQ:
          START_TIMING(steal_wait);
          MPI_Wait(&rq, &status);
          END_TIMING(steal_wait, elapsed_times[COMM_GLOBAL_TIME]);
          inter_node_work_steal_victim(availableWork, st, my_color, tok, comms);
          break;
        }
      }
      START_TIMING(win_sync);
      MPI_Win_sync(term_win);
      END_TIMING(win_sync, elapsed_times[COMM_GLOBAL_TIME]);
    }
    START_TIMING(send_wait);
    MPI_Wait(&rq, &status);
    END_TIMING(send_wait, elapsed_times[COMM_GLOBAL_TIME]);

    if (length(loot) > 0) {
      availableWork = assign_and_update_availableWork(loot, num_workers,
                                                      combs_dispensed, comms);
      lootReceived = true;
    }
  }
}

static void send_poison_pill(int num_workers, const CommsStruct &comms) {
  WorkChunk poison{0, -2};
  for (int w = 1; w <= num_workers; ++w) {
    START_TIMING(poison);
    MPI_Send(&poison, sizeof(poison), MPI_BYTE, w, TAG_ASSIGN_WORK,
             comms.local_comm);
    END_TIMING(poison, elapsed_times[COMM_LOCAL_TIME]);
  }
}

static void node_leader_hierarchical(WorkChunk availableWork, int num_workers,
                                     LAMBDA_TYPE num_Comb,
                                     const CommsStruct &comms) {
  bool *global_done;
  MPI_Win term_win;
  MPI_Win_allocate(sizeof(bool), sizeof(bool), MPI_INFO_NULL, comms.global_comm,
                   &global_done, &term_win);
  *global_done = false;
  MPI_Win_lock_all(0, term_win);

  const int next_leader = (comms.global_rank + 1) % comms.num_nodes;
  Token tok = {WHITE, false};
  bool have_token = (comms.global_rank == 0);
  int my_color = WHITE;
  std::size_t leader_iter = 0;
  LAMBDA_TYPE combs_dispensed = CHUNK_SIZE * num_workers;
  while (true) {
    MPI_Win_sync(term_win);
    if (*global_done) {
      break;
    }
    MPI_Status st;
    int flag;
    // local probe
    START_TIMING(local_probe);
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comms.local_comm, &flag, &st);
    END_TIMING(local_probe, elapsed_times[COMM_LOCAL_TIME]);
    if (flag) {
      int tag = st.MPI_TAG;
      switch (tag) {
      case TAG_REQUEST_WORK:
        handle_local_work_steal(availableWork, st, combs_dispensed, comms);
        break;
      }
    }

    // global probe
    flag = 0;
    START_TIMING(global_probe);
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comms.global_comm, &flag, &st);
    END_TIMING(global_probe, elapsed_times[COMM_GLOBAL_TIME]);
    if (flag) {
      int tag = st.MPI_TAG;
      switch (tag) {
      case TAG_NODE_STEAL_REQ:
        inter_node_work_steal_victim(availableWork, st, my_color, tok, comms);
        break;
      case TAG_TOKEN:
        receive_token(tok, st, have_token, comms);
        break;
      }
    }
    try_forward_token_if_idle(availableWork, have_token, my_color, tok,
                              next_leader, term_win, comms);
    // If leader node is idle, initiate a steal request
    if (length(availableWork) <= 0 && !(*global_done)) {
      inter_node_work_steal_initiate(availableWork, st, num_workers, my_color,
                                     tok, term_win, global_done,
                                     combs_dispensed, comms);
    }

#ifdef ENABLE_PROFILE
    if ((leader_iter % PRINT_FREQ == 0)) {
      START_TIMING(print_leader);

      double now = MPI_Wtime();
      double total_outer_elapsed = now - gprog.dist_start_ts;
      double avg_outer_time =
          (gprog.dist_iters_completed > 0)
              ? (gprog.outer_time_sum / (double)gprog.dist_iters_completed)
              : 0.0;
      double inner_elapsed = now - gprog.inner_start_ts;

      std::size_t iter_display = gprog.dist_iters_completed + 1;

      LAMBDA_TYPE approx_per_node = num_Comb / comms.num_nodes;
      double inner_pct =
          100.0 * (double)combs_dispensed / (double)approx_per_node;

      /**
      printf("iter: %zu | cover: %lld/%lld | time: %.0f sec | avg_outer_time: "
             "%.0f sec "
             "||| inner_progress (combs dispensed): ~%lld/%lld (~%.0f%%) | "
             "tasks unclaimed [start -> end]: %lld -> %lld | "
             "inner_time: %.0f sec\n",
             iter_display, gprog.cover_count, gprog.total_tumor,
             total_outer_elapsed, avg_outer_time, (long long)combs_dispensed,
             (long long)approx_per_node, inner_pct, availableWork.start,
             availableWork.end, inner_elapsed); **/

      double comm_local_iter = iter_delta(COMM_LOCAL_TIME);
      double comm_global_iter = iter_delta(COMM_GLOBAL_TIME);
      double worker_idle_iter = iter_delta(WORKER_IDLE_TIME);
      double worker_time_iter = iter_delta(WORKER_TIME);
      double worker_run_direct_iter = iter_delta(WORKER_RUNNING_TIME);
      double worker_run_iter =
          (worker_run_direct_iter > 0.0)
              ? worker_run_direct_iter
              : std::max(0.0, worker_time_iter - worker_idle_iter);

      printf("iter: %zu | cover: %lld/%lld | time: %.0f sec | avg_outer_time: "
             "%.0f sec "
             "||| inner_progress (combs dispensed): ~%lld/%lld (~%.0f%%) | "
             "tasks unclaimed [start -> end]: %lld -> %lld | "
             "inner_time: %.0f sec | "
             "comm_local(iter): %.0f sec | comm_global(iter): %.0f sec | "
             "worker_idle(iter): %.0f sec | worker_run(iter): %.0f sec\n",
             iter_display, gprog.cover_count, gprog.total_tumor,
             total_outer_elapsed, avg_outer_time, (long long)combs_dispensed,
             (long long)approx_per_node, inner_pct, availableWork.start,
             availableWork.end, inner_elapsed, comm_local_iter,
             comm_global_iter, worker_idle_iter, worker_run_iter);
      fflush(stdout);

      END_TIMING(print_leader, elapsed_times[EXCLUDE_TIME]);
    }
    ++leader_iter;
#endif
  }
  // Poison the workers
  send_poison_pill(num_workers, comms);
  START_TIMING(free_window);
  MPI_Win_unlock_all(term_win);
  MPI_Win_free(&term_win);
  END_TIMING(free_window, elapsed_times[COMM_GLOBAL_TIME]);
}

static void worker_hierarchical(int worker_local_rank, WorkChunk &myChunk,
                                double &localBestMaxF, int localComb[],
                                sets_t dataTable, SET *buffers,
                                CommsStruct &comms) {
  MPI_Status status;

  while (true) {
    process_lambda_interval(myChunk.start, myChunk.end, localComb,
                            localBestMaxF, dataTable, buffers, comms);
    START_TIMING(idle);
    char dummy;
    START_TIMING(local_steal);
    MPI_Send(&dummy, 1, MPI_BYTE, 0, TAG_REQUEST_WORK, comms.local_comm);
    MPI_Recv(&myChunk, sizeof(WorkChunk), MPI_BYTE, 0, TAG_ASSIGN_WORK,
             comms.local_comm, &status);
    END_TIMING(local_steal, elapsed_times[COMM_LOCAL_TIME]);

    if (length(myChunk) < 0) {
      END_TIMING(idle, elapsed_times[WORKER_IDLE_TIME]);
      break;
    }
    END_TIMING(idle, elapsed_times[WORKER_IDLE_TIME]);
  }
}

static inline void execute_hierarchical(int rank, int size_minus_one,
                                        LAMBDA_TYPE num_Comb,
                                        double &localBestMaxF, int localComb[],
                                        sets_t dataTable, SET *buffers,
                                        CommsStruct &comms) {

#ifdef ENABLE_PROFILE
  snapshot_iter_timers();
#endif

#ifdef ENABLE_PROFILE
  if (comms.global_rank == 0)
    gprog.inner_start_ts = MPI_Wtime();
#endif

  WorkChunk leaderRange = calculate_node_range(num_Comb, comms);
  const int num_workers = comms.local_size - 1;

  if (comms.is_leader) {
    leaderRange.start += (CHUNK_SIZE * num_workers);
    START_TIMING(leader_time);
    node_leader_hierarchical(leaderRange, num_workers, num_Comb, comms);
    END_TIMING(leader_time, elapsed_times[MASTER_TIME]);
  } else {
    const int worker_id = comms.local_rank - 1;
    WorkChunk myChunk =
        calculate_worker_range(leaderRange, worker_id, num_workers);

    START_TIMING(worker_time);
    worker_hierarchical(comms.local_rank, myChunk, localBestMaxF, localComb,
                        dataTable, buffers, comms);
    END_TIMING(worker_time, elapsed_times[WORKER_TIME]);
  }
}

#else // Not using hierachical Allreduce

#define ALL_REDUCE_FUNC MPI_Allreduce
#define EXECUTE execute_role

#endif // ALL_REDUCE_HIERARCHICAL

//////////////////////////////  End Allreduce_hierarchical
/////////////////////////

static inline LAMBDA_TYPE calculate_initial_index(int num_workers) {
  return static_cast<LAMBDA_TYPE>(num_workers) * CHUNK_SIZE;
}

static inline void distribute_work(int num_workers, LAMBDA_TYPE num_Comb,
                                   LAMBDA_TYPE &next_idx) {
  while (next_idx < num_Comb) {
    MPI_Status status;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

    if (flag) {
      int workerRank = status.MPI_SOURCE;
      char message;
      MPI_Recv(&message, 1, MPI_CHAR, workerRank, 1, MPI_COMM_WORLD, &status);

      if (message == 'a') {
        MPI_Send(&next_idx, 1, MPI_LONG_LONG_INT, workerRank, 2,
                 MPI_COMM_WORLD);
        std::cout << "[At time epoch: " << MPI_Wtime() << "s] "
                  << "Sent combos [" << next_idx << "–"
                  << (std::min(next_idx + CHUNK_SIZE, num_Comb) - 1)
                  << "] to rank " << workerRank << ". "
                  << (num_Comb - std::min(next_idx + CHUNK_SIZE, num_Comb))
                  << " combos left." << std::endl;
        next_idx += CHUNK_SIZE;
      }
    }
  }
}

static inline void master_process(int num_workers, LAMBDA_TYPE num_Comb) {
  LAMBDA_TYPE next_idx = calculate_initial_index(num_workers);
  distribute_work(num_workers, num_Comb, next_idx);

  LAMBDA_TYPE termination_signal = -1;
  for (int workerRank = 1; workerRank <= num_workers; ++workerRank) {
    MPI_Send(&termination_signal, 1, MPI_LONG_LONG_INT, workerRank, 2,
             MPI_COMM_WORLD);
  }
}

static inline LAMBDA_TYPE nCr(int n, int r) {
  if (r > n)
    return 0;
  if (r == 0 || r == n)
    return 1;
  if (r > n - r)
    r = n - r; // Because C(n, r) = C(n, n-r)

  LAMBDA_TYPE result = 1;
  for (int i = 1; i <= r; ++i) {
    result *= (n - r + i);
    result /= i;
  }
  return result;
}

static void outputFileWriteError(std::ofstream &outfile) {

  if (!outfile) {
    std::cerr << "Error: Could not open output file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

std::pair<LAMBDA_TYPE, LAMBDA_TYPE>
calculate_initial_chunk(int rank, LAMBDA_TYPE num_Comb,
                        LAMBDA_TYPE chunk_size) {
  LAMBDA_TYPE begin = (rank - 1) * chunk_size;
  LAMBDA_TYPE end = std::min(begin + chunk_size, num_Comb);
  return {begin, end};
}

static inline LambdaComputed compute_lambda_variables(LAMBDA_TYPE lambda,
                                                      int totalGenes) {
  LambdaComputed computed;
  computed.j = static_cast<int>(std::floor(std::sqrt(0.25 + 2 * lambda) + 0.5));
  if (computed.j > totalGenes - (NUMHITS - 2)) {
    computed.j = -1;
    return computed;
  }
  computed.i = lambda - (computed.j * (computed.j - 1)) / 2;
  return computed;
}

static void write_output(int rank, std::ofstream &outfile,
                         const int globalBestComb[], double F_max,
                         const LAMBDA_TYPE boundCounts[NUMHITS],
                         LAMBDA_TYPE totalCombPossible) {
  outfile << "(";
  for (size_t idx = 0; idx < NUMHITS; ++idx) {
    outfile << globalBestComb[idx];
    if (idx != NUMHITS - 1) {
      outfile << ", ";
    }
  }
#ifdef ENABLE_PROFILE
  outfile << ")  F-max = " << F_max << ", Prune distribution (level 0 ... "
          << NUMHITS - 1 << "): ";
  for (int i = 0; i < NUMHITS; ++i) {
    outfile << boundCounts[i]
            << (i == NUMHITS - 1 ? " <-- Last index is number Unpruned" : ", ");
  }
  outfile << " out of " << totalCombPossible << " combinations." << std::endl;
#else
  outfile << ")  F-max = " << F_max << std::endl;
#endif // ENABLE_PROFILE
}

static void max_f_with_comb(void *in, void *inout, int *len,
                            MPI_Datatype *type) {
  const MPIResultWithComb *in_vals = static_cast<const MPIResultWithComb *>(in);
  MPIResultWithComb *inout_vals = static_cast<MPIResultWithComb *>(inout);

  for (int i = 0; i < *len; i++) {
    if (in_vals[i].f > inout_vals[i].f) {
      inout_vals[i].f = in_vals[i].f;
      for (int j = 0; j < NUMHITS; j++) {
        inout_vals[i].comb[j] = in_vals[i].comb[j];
      }
    }
  }
}

static MPI_Op create_max_f_with_comb_op(MPI_Datatype MPI_RESULT_WITH_COMB) {
  MPI_Op MPI_MAX_F_WITH_COMB;
  MPI_Op_create(&max_f_with_comb, 1, &MPI_MAX_F_WITH_COMB);
  return MPI_MAX_F_WITH_COMB;
}

static MPI_Datatype create_mpi_result_with_comb_type() {
  MPI_Datatype MPI_RESULT_WITH_COMB;

  const int nitems = 2;
  int blocklengths[] = {1, NUMHITS};
  MPI_Datatype types[] = {MPI_DOUBLE, MPI_INT};

  MPI_Aint offsets[nitems];
  offsets[0] = offsetof(MPIResultWithComb, f);
  offsets[1] = offsetof(MPIResultWithComb, comb);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types,
                         &MPI_RESULT_WITH_COMB);
  MPI_Type_commit(&MPI_RESULT_WITH_COMB);

  return MPI_RESULT_WITH_COMB;
}

static inline void process_lambda_interval(LAMBDA_TYPE startComb,
                                           LAMBDA_TYPE endComb,
                                           int bestCombination[], double &maxF,
                                           sets_t &dataTable, SET *buffers,
                                           CommsStruct &comms) {

  const int totalGenes = dataTable.numRows;
  const double alpha = 0.1;
  int localComb[NUMHITS] = {0};

  for (LAMBDA_TYPE lambda = startComb; lambda <= endComb; ++lambda) {

    LambdaComputed computed = compute_lambda_variables(lambda, totalGenes);
    if (computed.j < 0) {
      continue;
    }

    SET rowI =
        GET_ROW(dataTable.tumorData, computed.i, dataTable.tumorRowUnits);
    SET rowJ =
        GET_ROW(dataTable.tumorData, computed.j, dataTable.tumorRowUnits);
    SET_INTERSECT(buffers[0], rowI, rowJ, dataTable.tumorRowUnits);

#ifdef BOUND
    if (SET_IS_EMPTY(buffers[0], dataTable.tumorRowUnits)) {
      INCREMENT_BOUND_LEVEL(1);
      continue;
    }
#endif

    localComb[0] = computed.i;
    localComb[1] = computed.j;

#if NUMHITS == 2
    const int TP = SET_COUNT(buffers[0], dataTable.tumorRowUnits);

    SET normalRows[2];
    normalRows[0] =
        GET_ROW(dataTable.normalData, localComb[0], dataTable.normalRowUnits);
    normalRows[1] =
        GET_ROW(dataTable.normalData, localComb[1], dataTable.normalRowUnits);
    SET_INTERSECT_N(buffers[0], normalRows, 2, dataTable.normalRowUnits);

    const int coveredNormal = SET_COUNT(buffers[0], dataTable.normalRowUnits);
    const int TN = (int)dataTable.numNormal - coveredNormal;
    const double F =
        (alpha * TP + TN) / (dataTable.numTumor + dataTable.numNormal);

    if (F >= maxF) {
      maxF = F;
      bestCombination[0] = localComb[0];
      bestCombination[1] = localComb[1];
    }
#ifdef BOUND
    INCREMENT_BOUND_LEVEL(NUMHITS - 1);
#endif
    continue;
#endif

    int indices[NUMHITS];
    indices[0] = computed.i;
    indices[1] = computed.j;
    indices[2] = computed.j + 1;
    int level = 2;

    while (level >= 2) {
      int maxStart = totalGenes - (NUMHITS - (level + 1));
      if (indices[level] >= maxStart) {
        --level;
        if (level >= 2) {
          ++indices[level];
        }
        continue;
      }

      SET rowK =
          GET_ROW(dataTable.tumorData, indices[level], dataTable.tumorRowUnits);
      SET_INTERSECT(buffers[level - 1], buffers[level - 2], rowK,
                    dataTable.tumorRowUnits);
#ifdef BOUND
      if (SET_IS_EMPTY(buffers[level - 1], dataTable.tumorRowUnits)) {
        INCREMENT_BOUND_LEVEL(level);
        ++indices[level];
        continue;
      }
#endif
      localComb[level] = indices[level];

      if (level == NUMHITS - 1) {
        int TP = SET_COUNT(buffers[NUMHITS - 2], dataTable.tumorRowUnits);
        SET normalRows[NUMHITS];
        for (int idx = 0; idx < NUMHITS; ++idx) {
          normalRows[idx] = GET_ROW(dataTable.normalData, localComb[idx],
                                    dataTable.normalRowUnits);
        }
        SET_INTERSECT_N(buffers[NUMHITS - 2], normalRows, NUMHITS,
                        dataTable.normalRowUnits);
        int coveredNormal =
            SET_COUNT(buffers[NUMHITS - 2], dataTable.normalRowUnits);
        int TN = (int)dataTable.numNormal - coveredNormal;
        double F =
            (alpha * TP + TN) / (dataTable.numTumor + dataTable.numNormal);
        if (F >= maxF) {
          maxF = F;
          for (int k = 0; k < NUMHITS; ++k)
            bestCombination[k] = localComb[k];
        }
        INCREMENT_BOUND_LEVEL(NUMHITS - 1);
        ++indices[level];
      } else {
        ++level;
        indices[level] = indices[level - 1] + 1;
      }
    }
  }
}

static bool process_and_communicate(int rank, LAMBDA_TYPE num_Comb,
                                    double &localBestMaxF, int localComb[],
                                    LAMBDA_TYPE &begin, LAMBDA_TYPE &end,
                                    MPI_Status &status, sets_t dataTable,
                                    SET *buffers, CommsStruct &comms) {
  START_TIMING(run_time);
  process_lambda_interval(begin, end, localComb, localBestMaxF, dataTable,
                          buffers, comms);
  END_TIMING(run_time, elapsed_times[WORKER_RUNNING_TIME]);
  char signal = 'a';
  MPI_Send(&signal, 1, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
  LAMBDA_TYPE next_idx;
  MPI_Recv(&next_idx, 1, MPI_LONG_LONG_INT, 0, 2, MPI_COMM_WORLD, &status);

  if (next_idx == -1) {
    return false;
  }

  begin = next_idx;
  end = std::min(begin + CHUNK_SIZE, num_Comb);
  return true;
}

static void worker_process(int rank, LAMBDA_TYPE num_Comb,
                           double &localBestMaxF, int localComb[],
                           sets_t dataTable, SET *buffers, CommsStruct &comms) {
  std::pair<LAMBDA_TYPE, LAMBDA_TYPE> chunk_indices =
      calculate_initial_chunk(rank, num_Comb, CHUNK_SIZE);

  LAMBDA_TYPE begin = chunk_indices.first;
  LAMBDA_TYPE end = chunk_indices.second;

  MPI_Status status;

  while (end <= num_Comb) {
    bool has_next =
        process_and_communicate(rank, num_Comb, localBestMaxF, localComb, begin,
                                end, status, dataTable, buffers, comms);
    if (!has_next) {
      break;
    }
  }
}

static void execute_role(int rank, int size_minus_one, LAMBDA_TYPE num_Comb,
                         double &localBestMaxF, int localComb[],
                         sets_t dataTable, SET *buffers, CommsStruct &comms) {
  if (rank == 0) {
#ifdef ENABLE_PROFILE
    gprog.inner_start_ts = MPI_Wtime();
#endif
    START_TIMING(master_proc);
    master_process(size_minus_one, num_Comb);
    END_TIMING(master_proc, elapsed_times[MASTER_TIME]);
  } else {
    START_TIMING(worker_proc);
    worker_process(rank, num_Comb, localBestMaxF, localComb, dataTable, buffers,
                   comms);
    END_TIMING(worker_proc, elapsed_times[WORKER_TIME]);
  }
}

static inline void initialize_local_comb_and_f(double &f, int localComb[]) {
  f = 0;
  for (int i = 0; i < NUMHITS; ++i) {
    localComb[i] = -1;
  }
}

static MPIResultWithComb create_mpi_result(double f, const int comb[]) {
  MPIResultWithComb result;
  result.f = f;
  for (int i = 0; i < NUMHITS; ++i) {
    result.comb[i] = comb[i];
  }
  return result;
}

static void extract_global_comb(int globalBestComb[],
                                const MPIResultWithComb &globalResult) {
  for (int i = 0; i < NUMHITS; ++i) {
    globalBestComb[i] = globalResult.comb[i];
  }
}

void distribute_tasks(int rank, int size, const char *outFilename,
                      sets_t dataTable, CommsStruct &comms) {
  int Nt = dataTable.numTumor;
  int numGenes = dataTable.numRows;

  size_t tumorBits = dataTable.numTumor;
  size_t tumorUnits = dataTable.tumorRowUnits;
  size_t maxUnits = std::max(dataTable.tumorRowUnits, dataTable.normalRowUnits);
  SET buffers[NUMHITS - 1];

  for (int i = 0; i < NUMHITS - 1; i++) {
    SET_NEW(buffers[i], maxUnits);
  }

  MPI_Datatype MPI_RESULT_WITH_COMB = create_mpi_result_with_comb_type();
  MPI_Op MPI_MAX_F_WITH_COMB = create_max_f_with_comb_op(MPI_RESULT_WITH_COMB);

  LAMBDA_TYPE num_Comb = nCr(numGenes, 2);
  LAMBDA_TYPE totalCombPossible = nCr(numGenes, NUMHITS);

  SET droppedSamples;
  SET_NEW(droppedSamples, tumorUnits);

  std::ofstream outfile;
  if (rank == 0) {
    outfile.open(outFilename);
    outputFileWriteError(outfile);
  }

#ifdef ENABLE_PROFILE
  gprog.dist_iters_completed = 0;
  gprog.cover_count = 0;
  gprog.total_tumor = (long long)dataTable.numTumor;
  gprog.dist_start_ts = MPI_Wtime();
  gprog.outer_time_sum = 0.0;
  gprog.inner_start_ts = 0.0;
#endif

  while (
      !CHECK_ALL_BITS_SET(droppedSamples, tumorBits, dataTable.tumorRowUnits)) {
#ifdef ENABLE_PROFILE
    double outer_iter_start = MPI_Wtime();
#endif

    std::fill(std::begin(bound_level_counts), std::end(bound_level_counts), 0);
    double localBestMaxF;
    int localComb[NUMHITS];
    initialize_local_comb_and_f(localBestMaxF, localComb);

    EXECUTE(rank, size - 1, num_Comb, localBestMaxF, localComb, dataTable,
            buffers, comms);

    MPIResultWithComb localResult = create_mpi_result(localBestMaxF, localComb);
    MPIResultWithComb globalResult = {};
    ALL_REDUCE_FUNC(&localResult, &globalResult, 1, MPI_RESULT_WITH_COMB,
                    MPI_MAX_F_WITH_COMB, comms);
    int globalBestComb[NUMHITS];
    extract_global_comb(globalBestComb, globalResult);

    START_TIMING(metrics_time);
    long long globalBoundCounts[NUMHITS] = {0};
    ALL_REDUCE_FUNC(bound_level_counts, globalBoundCounts, NUMHITS,
                    MPI_LONG_LONG, MPI_SUM, comms);
    END_TIMING(metrics_time, elapsed_times[EXCLUDE_TIME]);

    SET intersectionSets[NUMHITS];
    for (int i = 0; i < NUMHITS; ++i) {
      intersectionSets[i] =
          GET_ROW(dataTable.tumorData, globalBestComb[i], tumorUnits);
    }

    SET_INTERSECT_N(buffers[NUMHITS - 2], intersectionSets, NUMHITS,
                    tumorUnits);

    SET_UNION(droppedSamples, droppedSamples, buffers[NUMHITS - 2],
              dataTable.tumorRowUnits);

    UPDATE_SET_COLLECTION(dataTable.tumorData, buffers[NUMHITS - 2],
                          dataTable.numRows, dataTable.tumorRowUnits);

    Nt -= SET_COUNT(buffers[NUMHITS - 2], dataTable.tumorRowUnits);

    if (rank == 0) {
      write_output(rank, outfile, globalBestComb, globalResult.f,
                   globalBoundCounts, totalCombPossible);
    }
#ifdef ENABLE_PROFILE
    gprog.cover_count = SET_COUNT(droppedSamples, dataTable.tumorRowUnits);
    gprog.outer_time_sum += (MPI_Wtime() - outer_iter_start);
    gprog.dist_iters_completed += 1;
#endif
  }

  if (rank == 0)
    outfile.close();

  for (int i = 0; i < NUMHITS - 1; i++) {
    SET_FREE(buffers[i]);
  }
  SET_FREE(droppedSamples);

  MPI_Op_free(&MPI_MAX_F_WITH_COMB);
  MPI_Type_free(&MPI_RESULT_WITH_COMB);
}
