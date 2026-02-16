# P-DFS-Multihit-WSC

A C++ implementation of a **Pruned Depth-First Search (P-DFS)** pipeline for performing **Weighted Set Cover (WSC)** to identify gene combinations that contribute towards carcinogenesis.

---

## Requirements
- CMake 
- A C++17-capable compiler (GCC/Clang)
- Ninja
- **MPI runtime** for `mpirun`

---

## Build

This project is typically built with **Ninja**. If the repo already has a configured `build/` directory (i.e., it contains a `build.ninja` file), you can build directly:

```bash
cd build
ninja
```

If your `build/` directory is not configured yet (no `build.ninja`), configure it once with CMake and then build:

```bash
mkdir -p build
cd build
cmake -G Ninja ..
ninja
```

## Run

After compiling, an executable (named `run`) is produced by the build.

A typical MPI launch looks like:

```bash
mpirun -np <NUM_PROCESSES> <PATH_TO_EXECUTABLE> <PATH_TO_DATA_FILE> <PATH_TO_METRICS_FILE> <PATH_TO_OUTPUT_FILE>
```

### Arguments (generic)

- `<NUM_PROCESSES>`: number of MPI ranks to launch (e.g., 1, 2, 4, 8, …)
- `<PATH_TO_EXECUTABLE>`: path to the compiled binary 
- `<PATH_TO_DATA_FILE>`: your input dataset file
- `<PATH_TO_METRICS_FILE>`: **a file that will be created/overwritten** to store run metrics
- `<PATH_TO_OUTPUT_FILE>`: **a file that will be created/overwritten** to store the solver output 

### Example

```bash
mpirun -np <NUM_PROCESSES> ./run <DATA_FILE> <METRICS_FILE> <OUTPUT_FILE>
```



