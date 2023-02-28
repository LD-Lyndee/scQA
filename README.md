# scQA
Clustering scRNA-seq data via qualitative and quantitative analysis

## Installation:

Installation Tested on a Linux sever.

### Source:

```bash
git clone https://github.com/LD-Lyndee/scQA.git
  cd scQA
```

### Make File:

```shell
./make
```

## Quick Start

### 1. Prepare Input Matrix

Pre-processing with R provided as pre.R

### 2. Run scQA

#### Matrix Format

The input file must be in the form of a matrix of cells * genes, without names of rows and columns.

#### Parameter

There are two input parameters：

- **filename** Input matrix file. 
- **binning** Define number of bins, default: 6. "o" for using the default parameter.

```bash
./main xin_matrix.txt o
```

### 3. Output

- **labels.txt** Clustering result.




