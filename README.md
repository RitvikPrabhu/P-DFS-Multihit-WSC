# P-DFS-Multihit-WSC

A C++ implementation of a **Pruned Depth-First Search (P-DFS)** pipeline for performing **Weighted Set Cover (WSC)** to identify gene combinations that contribute towards carcinogenesis.

## Build (CMake)

### Requirements
- CMake 
- A C++17-capable compiler (GCC/Clang)
- Linux/macOS recommended

- ### Compile
```bash
# From repo root
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```
