# GLaMST
This algorithm is developed to reconstruct lineage tree from subset of tree nodes.

# Download 
GLaMST can be downloaded from https://github.com/xysheep/GLaMST/releases

# Usage of Source code
A simple demonstration of how to use the source code is avaliable at https://xysheep.github.io/GLaMST . If Matlab is not avaliable to you, you can use the pre-compiled version. 

# Usage of pre-compiled excutible
To use this executable, the MATLAB Runtime 9.2 is required. It can be downloaded from [here](https://www.mathworks.com/products/compiler/matlab-runtime.html)
```
GLaMST.exe fastafilename [...options]
```
# Arguments
## Required Arguments
- fastafilename: name of one fasta file. The first sequence in the fasta file should be root (G. L.). **Sequences in this fasta file should be original assembled sequences**, not multiple alignment results. 
## Optional Arguments
- rewire: 1 means do rewire while 0 means not do rewire [default 1]
