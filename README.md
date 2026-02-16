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

## Output format

The solver writes an output file containing **one candidate combination per line**, formatted as a parenthesized, comma-separated tuple of **gene indices**:

```
(i1, i2, i3, ..., ik)
(i1, i2, i3, ..., ik)
...
```

- The **tuple length = k**, where **k is the number of hits** you ran (e.g., `k=5` → 5 numbers per line).
- Each number is an **index pointer** into the gene list derived from your input data preprocessing (i.e., the file stores indices, not gene names).
- Interpreting a line: it represents the **set of k genes** (by index) forming one k-hit combination.

---

## Convert indices to gene names

To translate index tuples into actual gene names, use:

- `utils/convertIndexToGeneName.py`

### Usage

```bash
python3 utils/convertIndexToGeneName.py <RESULT_FILE> <GENE_MAP_FILE>
```

- `<RESULT_FILE>`: the solver output file containing tuples like `(i1, i2, ..., ik)`
- `<GENE_MAP_FILE>`: the gene map produced during data preprocessing/merging (when `utils/process_gene_data.py` was run) that defines the mapping from **index → gene name**

This produces a human-readable view of the results with gene names instead of indices (exact formatting depends on the script).

---

## Verify tumor coverage

In theory, taking the **union** of all k-hit sets in the result file should allow you to cover all tumor patients (subject to preprocessing, run parameters, and dataset constraints). You can verify this using:

- `utils/verifyAccuracy.py`

### Usage

```bash
python3 utils/verifyAccuracy.py <MATRIX_FILE> <GENE_SAMPLE_LIST_FILE> <RESULT_FILE>
```

- `<MATRIX_FILE>`: the preprocessed mutation matrix produced during raw-data preprocessing (**before** the merge step)
- `<GENE_SAMPLE_LIST_FILE>`: the gene-sample list produced during preprocessing (**before** the merge step)
- `<RESULT_FILE>`: the solver output file containing `(i1, i2, ..., ik)` tuples

The script reports whether the result file achieves the expected tumor coverage and typically identifies any uncovered samples (exact reporting depends on the script).

---




