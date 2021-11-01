/*
 *  CRSmatrix.h
 *  CRS matrix functions and Jacobi method
 *
 */

#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

/* sparse matrix structure using row-compressed storage
        Contents:
                value: sparse matrix values
                collnd: column index associated with each matrix value
                rowPtr: keeps count of the number of values in each row and the number of rows
                rowNum: number of rows
                colNum: number of columns
*/
class CRSmatrix {
protected:
        vector<double> value;
        vector<int> rowPtr;
        vector<int> colInd;
        int rowNum;
        int colNum;
public:
        //CRSmatrix constructor with value, row point, and column index vectors
        //Indices must be 0 based
        CRSmatrix(vector<double> v, vector<int> r, vector<int> c);

        //CRSmatrix constructor with local file paths containing the data
        //Indices must be 0 based
        CRSmatrix(string valueAddress, string rowPtrAddress, string colIndAddress);

        //CRSmatrix constructor with a mtx file
        //Indices must be 1 based (mtx files should be 1 based)
        CRSmatrix(string mtxFile);

        //CRSmatrix constructor with row number of column number, values will be filled with 0
        CRSmatrix(int r, int c);

        //get row and col size
        int getRowSize ();
        int getColSize ();

        //return matrix value at [i][j]
        double retrieveElement (int i, int j);

        //insert value x at position [i][j]
        void changeValue (double x, int i, int j);

        //multiply matrix with a vector x and return the product
        vector<double> productAx(vector<double> x);

        //delete value at position [i][j]
        void deleteValue(int i, int j);

        //print CRSmatrix
        void printA();
} ;

/* finds the vector norm
        Parameters:
                x: any vector
        Return:
                ||x||2
*/
double vectorNorm(vector<double> x);

/* Jacobi CRS matrix solver for Ax = b
        Parameters:
                A: a CRS matrix
                b: a vector
        Return:
                the solution to Ax = b
*/
vector<double> Jacobi(CRSmatrix A, vector<double> b);

#endif //CRSMATRIX_H