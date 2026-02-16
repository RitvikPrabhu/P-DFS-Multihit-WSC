# P-DFS-Multihit-WSC

A C++ implementation of a **Pruned Depth-First Search (P-DFS)** pipeline for performing **Weighted Set Cover (WSC)** to identify gene combinations that contribute towards carcinogenesis.

---

## Requirements
- CMake 
- A C++17-capable compiler (GCC/Clang)
- Ninja
- **MPI runtime** for `mpirun`
- Python 3.10 or later

---

## Data source (TCGA)

This codebase is designed to work with mutation data from **The Cancer Genome Atlas (TCGA)**:
- TCGA overview: https://www.cancer.gov/ccg/research/genome-sequencing/tcga
- Practical access is typically via the **Genomic Data Commons (GDC)** using the **GDC Data Portal** + **GDC client**.

### Download mutations (GDC client workflow)

1. In the **GDC Data Portal**, select the TCGA project(s) you want and add the corresponding mutation files (MAF) to your cart.
2. Download the **manifest** from the cart.
3. Use the **GDC client** to download the files:

```bash
gdc-client download -m <MANIFEST_FILE> -d <OUTPUT_DIR>
```

A common organization is one directory per cancer type / TCGA project:

```
<OUTPUT_DIR>/
  <TCGA_PROJECT_1>/   # MAF files
  <TCGA_PROJECT_2>/   # MAF files
  ...
```

---
## Data processing

The repository includes a preprocessing script:

- `utils/preprocessing_maf.py`

This script:
- reads all MAF files under a GDC download directory
- filters to functionally relevant variants
- separates tumor vs matched-normal calls
- emits intermediate files used by downstream merge/packaging steps

### Stage 1: Run MAF preprocessing

```bash
python3 utils/preprocessing_maf.py <GDC_DOWNLOAD_DIR>
```

**Input**
- `<GDC_DOWNLOAD_DIR>`: directory containing the GDC-downloaded MAF subfolders (the script scans `"<GDC_DOWNLOAD_DIR>/*/*.maf.*"`)

**Variant filtering**
- Keeps: `Missense_Mutation`, `Nonsense_Mutation`, `Frame_Shift_Ins`, `Frame_Shift_Del`,
  `In_Frame_Ins`, `In_Frame_Del`, `Splice_Site`, `Translation_Start_Site`, `Nonstop_Mutation`
- Normalizes sample barcodes to the first 12 characters.
- Treats a call as mutated if the tumor/normal allele differs from the reference allele.

### 2) Outputs produced by `preprocessing_maf.py`

The script writes **two files** (to the current working directory):

1. **Tumor mutation matrix file**
   - Name pattern:
     - `Tumor_matrix_attempt_2_<DATASET_TAG>.txt`
   - Contents:
     - a tab-separated table with gene/sample indices and a binary mutation indicator.
     - includes `Gene` and `Sample` (barcode) columns for traceability.

2. **Normal mutation list file**
   - Name pattern:
     - `Normal_list_attempt_2_<DATASET_TAG>.txt`
   - Contents:
     - a tab-separated two-column list containing `(Gene, SampleIndex)` pairs for normal calls.
     - sample indices are aligned to tumor indices when possible. Unmatched normals are assigned new indices.

> `<DATASET_TAG>` is derived from the input directory name inside the script.

---

### Stage 2: Create the combined solver input (bitwise + sorted)

After Stage 1, run:

- `utils/process_gene_data.py`

This script:
- combines the processed tumor/normal artifacts into a **bitwise representation**
- **sorts genes from most sparse → least sparse**
- writes the final **combined dataset file** used by `./run`
- produces a **gene map** used later to convert output indices back to gene names

#### Usage

```bash
python3 utils/process_gene_data.py <TUMOR_MATRIX_FILE> <GENE_SAMPLE_LIST_FILE> <COMBINED_OUTPUT_FILE>
```

- `<TUMOR_MATRIX_FILE>`: tumor mutation matrix produced during preprocessing
- `<GENE_SAMPLE_LIST_FILE>`: gene-sample list produced during preprocessing (**before** merge)
- `<COMBINED_OUTPUT_FILE>`: output file to be created (this is what you pass as `<PATH_TO_DATA_FILE>` to `./run`)

Additional output generated:
- `<COMBINED_OUTPUT_FILE>.gene_map` (index → gene name mapping)

#### What the input matrix looks like (conceptual)

Rows are genes; columns are samples (tumor + normal). Entries are binary mutation indicators.

```
        TS1  TS2  TS3  NS1  NS2
G1       1    0    1    0    0
G2       0    1    0    0    1
G3       1    1    0    0    0
```

(Here `TS*` are tumor samples and `NS*` are normal samples. The real files will use TCGA-style sample IDs.)

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

---

## Configuring number of hits and bounded vs exhaustive search

Two key settings are compile-time configuration options controlled by CMake cache variables:

1. **Number of hits (`NUMHITS`)**  
   In `CMakeLists.txt`, this is defined as:
   - `set(NUMHITS "4" CACHE STRING "Number of hits in the pattern (e.g. 3, 4, 5, 6...)")`

2. **Bounded/pruned vs exhaustive (`BOUND`)**  
   In `CMakeLists.txt`, this is defined as:
   - `option(BOUND "Enable DFS" ON)`

### Easiest way (recommended): set CMake variables on configure

Instead of editing `CMakeLists.txt`, pass these values when you run `cmake`:

```bash
mkdir -p build
cd build
cmake -G Ninja .. -DNUMHITS=<K> -DBOUND=ON
ninja
```

To switch to exhaustive mode, set `BOUND=OFF`:

```bash
cmake -G Ninja .. -DNUMHITS=<K> -DBOUND=OFF
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

- `<MATRIX_FILE>`: the tumor mutation matrix produced during preprocessing (**before** the merge step)
- `<GENE_SAMPLE_LIST_FILE>`: gene-sample list produced during preprocessing (**before** merge)
- `<RESULT_FILE>`: the solver output file containing `(i1, i2, ..., ik)` tuples

The script reports whether the result file achieves the expected tumor coverage and typically identifies any uncovered samples (exact reporting depends on the script).

---








