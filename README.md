# Sparse-Matrix
c++ compressed row storage sparse matrix ver 2.5

The CRS (Compressed Row Storage) Sparse Matrix was originally a programing assignment for ECE 4960. All assignments in this course were written from scratch. In the end of the course, I made major modifications in CRSmatrix version 2.0 and used it in my final project (Spline Fitting, which can be found here https://github.com/fanyazhi/Spline-Fitting). 

In version 2.5, I added in examples, made the style better, and added more documentation.

## Features:
* Sparse matrix with row-compressed storage: initializations, retrieve element, insert or change element, multiplication of sparse matrix with vector, print matrix
* Find vector norm
* Solve for Ax = b using the Jacobian Iterative Method

## A little bit of the math background:

**CRS Sparse Matrix**

"Sparse matrices" are matrices with most of their values zero. A of practical problems involve sparse matrices. Those problems require handling large sparse matrices efficiently in terms of memory usage and computation. A sparse matrix can be represented by an array of row index, an array of col index, and an array of the values. 

A further improvement on the sparse matrix is the CRS (Compressed Row Storage) Sparse Matrix. In the CRS representation, the row index is omitted, instead we use the turning values as the row pointers. This further reduce the amount of memory needed to store the matrix.

For example: <br />
|1 2 0| <br />
|3 0 0| <br />
|4 5 6| <br />
can be represented by <br />
value{1, 2, 3, 4, 5, 6} <br />
rowPtr{0, 2, 3, 6} <br />
colInd{0, 1, 0, 0, 1, 2} <br />

**Jacobian Iterative Solver**

When we are solving a large matrix A*x = b, direct methods involving fill-ins are too expensive. We instead need to use iterative solvers, such as the Jacobian Iterative Solver. Iterative solvers provide an approximation of the solution after a finite number of steps.

