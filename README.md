Parallel Direct N-Body
=====

This is a dense kernel matrix algorithms library which aids in the research, development, and use of parallel methods for systems of equations of the form:

![equation](http://latex.codecogs.com/gif.latex?r_i%3D%5Csum_jK%28t_i%2Cs_j%29%5C%2Cc_j)<br/>
where<br/>
![equation](http://latex.codecogs.com/gif.latex?K) is the _kernel_ generating the elements of the matrix,<br/>
![equation](http://latex.codecogs.com/gif.latex?s_j) are the _sources_ of the kernel,<br/>
![equation](http://latex.codecogs.com/gif.latex?c_j) are the _charges_ of the sources,<br/>
![equation](http://latex.codecogs.com/gif.latex?t_i) are the _targets_ of the kernel (which may be equivalent to the sources),<br/>
![equation](http://latex.codecogs.com/gif.latex?r_i) are the _results_.<br/>

This is a kernel-matrix equation. Matrices of this form can be found in a wide variety of fields include physics, statistics, and machine learning. This library provides an STL-like interface for direct evaluation of dense kernel matrix-vector products using communication optimal parallel algorithms.

Primary Authors:
* Wesley Chen
* Cris Cecka (ccecka@seas.harvard.edu)

Dependencies:
* C++ compiler with C++11 support.
* MPI

Building:
* 'make'

Supported Flags:
* P2P_DECAY_ITERATOR={0.1}
** Find and decay contiguous iterators to pointers to exploit blocking and SMP.
* P2P_BLOCK_SIZE=###
** Maximum block size of the recursive P2P blocked evaluation. (Deprecate?)
* P2P_NUM_THREADS=###
** Number of SMB threads to use in the recursive P2P blocked evaluation.
