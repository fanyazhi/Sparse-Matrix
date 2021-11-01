/*
 *  CRSmatrix.h
 *  CRS matrix functions and Jacobi method
 *
 */

#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

/* sparse matrix structure using row-compressed storage
        Contents:
                value: sparse matrix values
                collnd: column index associated with each matrix value
                rowPtr: keeps count of the number of values in each row and the number of rows
                rowNum: number of rows
                colNum: number of columns
*/
class CRSmatrix {
   private:
    int rowNum;
    int colNum;

   public:
    std::vector<double> value;
    std::vector<int> rowPtr;
    std::vector<int> colInd;

    //CRSmatrix constructor with value, row point, and column index vectors
    //Indices must be 0 based
    CRSmatrix(std::vector<double> v, std::vector<int> r, std::vector<int> c);

    //CRSmatrix constructor with local file paths containing the data
    //Indices must be 0 based
    CRSmatrix(std::string valueAddress, std::string rowPtrAddress, std::string colIndAddress);

    //CRSmatrix constructor with row number of column number, values will be filled with 0
    CRSmatrix(int r, int c);

    //get row and col size
    int getRowSize();
    int getColSize();

    //return matrix value at [i][j]
    double retrieveElement(unsigned int i, unsigned int j);

    //insert value x at position [i][j]
    void changeValue(double x, unsigned int i, unsigned int j);

    //multiply matrix with a vector x and return the product
    std::vector<double> productAx(std::vector<double> x);

    //delete value at position [i][j]
    void deleteValue(unsigned int i, unsigned int j);

    //print CRSmatrix
    void printA();
};

/* finds the vector norm
        Parameters:
                x: any vector
        Return:
                ||x||2
*/
double vectorNorm(std::vector<double> x);

/* Jacobi CRS matrix solver for Ax = b
        Parameters:
                A: a CRS matrix
                b: a vector
        Return:
                the solution to Ax = b
*/
std::vector<double> Jacobi(CRSmatrix A, std::vector<double> b);

#endif  //CRSMATRIX_H