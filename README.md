# P-DFS-Multihit-WSC

A C++ implementation of a **Pruned Depth-First Search (P-DFS)** pipeline for performing **Weighted Set Cover (WSC)** to identify gene combinations that contribute towards carcinogenesis.

## Build 

### Requirements
- CMake 
- A C++17-capable compiler (GCC/Clang)
- Ninja

### Compile

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


