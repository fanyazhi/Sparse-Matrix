/*
 *  The example in memplus.mtx is a large sparse matrix downloaded from the 
 *  matrix market (http://math.nist.gov/MatrixMarket/extreme.html) called memplus (rank: 17758; 
 *  126,150 nonzero entries) generated from the SPICE simulation of memory circuits. 
 *  Note that the mtx file format index is 1 based. 
 *  
 */
#include <iostream>
#include <cmath>
#include <vector>
#include "CRSmatrix.h"

using namespace std;

int main() {
	/*
	 *  1. construct sparse matrix from file 
	 */
    CRSmatrix A("../examples/value.txt", "../examples/rowPtr.txt", "../examples/colInd.txt");
	CRSmatrix B("../examples/memplus.mtx");

	/*
	 *  2. initialize a vector b
	 *  in this example, b=(1.0, 0, 0, …,0)^T, T is the number of columns in the matrix
	 */
    vector<double> b (A.colNum);
    for (int i = 0; i < b.size(); i++) {
		b[i] = 0;
	}
    b[0] = 1;

	/*
	 *  Optional: print out A and b
	 */

	cout << "b = [ ";
	for (int i = 0; i < b.size(); i++) {
		cout << b[i] << ", ";
	}
	cout << " ] " << endl;

	A.printA();

	/*
	 *  3. solver for Ax = b with the Jacobi method
	 */

	cout << "calculating Ax = b..." << endl;
    vector<double> x (b.size());
    x = Jacobi(A, b);

	/*
	 *  Optional: print out the resulting x
	 */
	cout << "x = [ ";
	for (int i = 0; i < x.size(); i++) {
		cout << x[i] << " ";
	}
	cout << " ] " << endl;
    
    return 0;
}