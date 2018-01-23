# GLaMST
This algorithm is developed to reconstruct lineage tree of hyper-mutation from subset of tree nodes.

Download GLaMST from [this link](https://github.com/xysheep/GLaMST/releases)

# Usage
To use this executable, the MATLAB Runtime 9.2 is required. It can be downloaded from [here](https://www.mathworks.com/products/compiler/matlab-runtime.html)
```
LineageTree.exe fastafilename [...options]
```
# Arguments
## Required Arguments
- fastafilename: name of one fasta file. The first sequence in the fasta file should be root (G. L.). 
## Optional Arguments
- rewire: 1 means do rewire while 0 means not do rewire [default 1]
