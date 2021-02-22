# analysis_template
Template of analysis made with AnalysisTree package.

# Requirement

* CMake 3.15
* ROOT v6-20-04 built with C++17
* AnalysisTree built with C++17

# Building

Clone this repository with
```
  git clone https://github.com/mam-mih-val/analysis_template
```
Create build directory
```
  cd analysis_template
  mkdir build
  cd build
```
Source root environment
```
  source /path/to/root/install/bin/thisroot.sh
``` 
Export AnalysisTree library
```
  export AnalysisTree_DIR=/path/to/AnalysisTree/install/lib/cmake/AnalysisTree/
```
Build the project
```
  cmake ..
  make -j
```

# Usage
To use the program run
```
  ./analyse path/to/file.list
```
Example of file list you can find in "lists" directory
