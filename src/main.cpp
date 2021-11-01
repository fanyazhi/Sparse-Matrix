/*
 * This files shows an example of solving Ax = b using the Jaconbian method
 *  
 */
#include <cmath>
#include <iostream>
#include <vector>

#include "CRSmatrix.h"

int main() {
    /*
	 *  1. construct sparse matrix 
	 *     there are three ways to initialize a matrix:
	 *       - initialize from value, row point, and column index vectors in CRS format (0 based)
	 *       - initialize from three files containing value, row point, and column index vectors in CRS format (0 based)
	 *       - initialize an empty matrix of size r x c, then use method CRSmatrix::changeValue to add in values
	 */
    std::vector<double> value{-4, 1, 1, 4, -4, 1, 1, -4, 1, 1, -4, 1, 1, 1, -4};
    std::vector<int> rowPtr{0, 3, 6, 9, 12, 15};
    std::vector<int> colInd{0, 1, 4, 0, 1, 2, 1, 2, 3, 2, 3, 4, 0, 3, 4};

    CRSmatrix matrixFromVectors(value, rowPtr, colInd);
    CRSmatrix matrixFrom3Files("../examples/value.txt", "../examples/rowPtr.txt", "../examples/colInd.txt");

    CRSmatrix A = matrixFrom3Files;

    /*
	 *  2. initialize a vector b
	 *  in this example, b=(1.0, 0, 0, …,0)^T, T is the number of columns in the matrix
	 */
    std::vector<double> b(A.getColSize());
    for (unsigned int i = 0; i < b.size(); i++) {
        b[i] = 0;
    }
    b[0] = 1;

    /*
	 *  Optional: print out A and b
	 */
    std::cout << "b = [ ";
    for (unsigned int i = 0; i < b.size(); i++) {
        std::cout << b[i] << ", ";
    }
    std::cout << " ] " << std::endl;

    std::cout << "A = " << std::endl;
    A.printA();

    /*
	 *  3. solve for Ax = b with the Jacobi method
	 */
    std::cout << "calculating Ax = b..." << std::endl;
    std::vector<double> x(b.size());
    x = Jacobi(A, b);

    /*
	 *  Optional: print out the resulting x
	 *
	 *  Sanity check: result for the matrixFromVectors and matrixFrom3Files should be:
     *  x = [ -0.37931 -0.408347 -0.116152 -0.0562609 -0.108892  ]
	 */
    std::cout << "x = [ ";
    for (unsigned int i = 0; i < x.size(); i++) {
        std::cout << x[i] << " ";
    }
    std::cout << " ] " << std::endl;

    return 0;
}