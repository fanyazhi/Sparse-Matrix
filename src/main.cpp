/*
 *  The example in memplus.mtx is a large sparse matrix downloaded from the 
 *  matrix market (http://math.nist.gov/MatrixMarket/extreme.html) called memplus (rank: 17758; 
 *  126,150 nonzero entries) generated from the SPICE simulation of memory circuits. 
 *  Note that the mtx file format index is 1 based. 
 * 
 * 
 *  
 */
#include <iostream>
#include <cmath>
#include <vector>
#include "CRSmatrix.h"

using namespace std;

int main() {
	/*
	 *  1. construct sparse matrix 
	 *     there are three ways to initialize a matrix:
	 *       - initialize from value, row point, and column index vectors in CRS format (0 based)
	 *       - initialize from three files containing value, row point, and column index vectors in CRS format (0 based)
	 *       - initialize from a mtx file in nomal sparse matrix format (1 based)
	 */
	vector<double> value {-4,1,1, 4,-4,1, 1,-4,1, 1,-4,1, 1,1,-4};
    vector<int> rowPtr {0, 3, 6, 9, 12, 15};
    vector<int> colInd {0,1,4, 0,1,2, 1,2,3, 2,3,4, 0,3,4};

    CRSmatrix matrixFromVectors(value, rowPtr, colInd);
	// CRSmatrix matrixFromMTX("../examples/memplus.mtx");
	CRSmatrix matrixFrom3Files("../examples/value.txt", "../examples/rowPtr.txt", "../examples/colInd.txt");

	CRSmatrix A = matrixFrom3Files;

	/*
	 *  2. initialize a vector b
	 *  in this example, b=(1.0, 0, 0, …,0)^T, T is the number of columns in the matrix
	 */
    vector<double> b (A.getColSize());
    for (unsigned int i = 0; i < b.size(); i++) {
		b[i] = 0;
	}
    b[0] = 1;

	/*
	 *  Optional: print out A and b, do not print if A is too large
	 */
	std::cout << "b = [ ";
	for (unsigned int i = 0; i < b.size(); i++) {
		cout << b[i] << ", ";
	}
	std::cout << " ] " << endl;

	std::cout << "A = " << endl;
	A.printA();

	/*
	 *  3. solve for Ax = b with the Jacobi method
	 */
	cout << "calculating Ax = b..." << endl;
    vector<double> x (b.size());
    x = Jacobi(A, b);

	/*
	 *  Optional: print out the resulting x
	 *
	 *  Sanity check: result for the matrixFromVectors and matrixFrom3Files  
     *  should be x = [ -0.37931 -0.408347 -0.116152 -0.0562609 -0.108892  ]
	 */
	cout << "x = [ ";
	for (unsigned int i = 0; i < x.size(); i++) {
		cout << x[i] << " ";
	}
	cout << " ] " << endl;

    return 0;
}