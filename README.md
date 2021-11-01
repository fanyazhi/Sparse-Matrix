# Sparse-Matrix
c++ compressed row storage sparse matrix ver 2.5

The CRS (Compressed Row Storage) Sparse Matrix was originally the first major programing assignment from ECE 4960, one of the most useful course I took at Cornell in 2018. In the end of the course, I made some major modifications to version 2.0 and used it in my final project (Spline Fitting, which is pretty cool). 

In version 2.5, the goal is to make the style better, the api easier to use, and add more documentation.

## Features:
* Sparse matrix structure with row-compressed storage: initialization of a sparse matrix, retrieve element, insert or change element, multiplication of sparse matrix with vector, print matrix
* Find vector norm
* Solve for Ax = b using the Jacobian Method

## A little bit of math backgrounds:
**CRS Sparse Matrix**
"Sparse matrices" are matrices with most of their values zero. A of practical problems iinvolve sparse matrices. Those problems require handling large sparse matrices efficiently in terms of memory usage and computation. A sparse matrix can be represented by an array of row index, an array of col index, and an array of the values. 

A further improvement on the sparse matrix is the CRS (Compressed Row Storage) Sparse Matrix. In the CRS representation, the row index is omitted, instead we use the "turning point values" as the row pointers. 


For example: <br />
|1 2 0| <br />
|3 0 0| <br />
|4 5 6| <br />
can be represented by <br />
value{1, 2, 3, 4, 5, 6} <br />
rowPtr{0, 2, 3, 6} <br />
colInd{0, 1, 0, 0, 1, 2} <br />

**Jacobian Iterative Solver**
When we are solving a large matrix A*x = b, direct methods involving fill-in and LU decomposition are too expensive. We need to use iterative solvers, such as the Jacobian Iterative Solver. Iterative solvers provide an approximation of the solution after a finite number of steps.

