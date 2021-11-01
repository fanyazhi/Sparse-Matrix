/*
 *  CRSmatrix.cpp
 *  CRS matrix functions and the Jacobi method
 *
 */

#include "CRSmatrix.h"

CRSmatrix::CRSmatrix(std::vector<double> v, std::vector<int> r, std::vector<int> c) {
    value = v;
    rowPtr = r;
    colInd = c;

    //find total number of rows
    rowNum = rowPtr.size() - 1;

    //find total number of columns
    colNum = 0;
    for (unsigned int i = 0; i < colInd.size(); i++) {
        if (colInd[i] > colNum) {
            colNum = colInd[i];
        }
    }
    colNum++;
}

CRSmatrix::CRSmatrix(std::string valueFilePath, std::string rowPtrFilePath, std::string colIndFilePath) {
    //Read in files
    std::ifstream rowFile, valueFile, colFile;

    rowFile.open(rowPtrFilePath);
    valueFile.open(valueFilePath);
    colFile.open(colIndFilePath);

    // Check if the file is open
    if (!rowFile || !valueFile || !colFile) {
        std::cout << "Unable to open file" << std::endl;
        exit(1);
    }

    // Declare a variable to load the contents from the file
    std::string line = "";

    // Loading the value vector
    std::string myValue;
    while (valueFile.good()) {
        getline(valueFile, myValue, '\n');
        //convert string value to double value
        double myDoubleValue = atof(myValue.c_str());
        value.push_back(myDoubleValue);
    }

    // Loading the row pointer vector
    while (rowFile >> line) {
        int rowNum = stoi(line);
        rowPtr.push_back(rowNum);
    }

    // Loading the column indice vector
    while (colFile >> line) {
        int colIndex = stoi(line);
        colInd.push_back(colIndex);
    }

    //find total number of rows
    rowNum = rowPtr.size() - 1;

    //find total number of columns
    colNum = 0;
    for (unsigned int i = 0; i < colInd.size(); i++) {
        if (colInd[i] > colNum) colNum = colInd[i];
    }
    colNum++;
}

CRSmatrix::CRSmatrix(int r, int c) {
    rowNum = r;
    colNum = c;
}

int CRSmatrix::getRowSize() {
    return rowNum;
}

int CRSmatrix::getColSize() {
    return colNum;
}

double CRSmatrix::retrieveElement(unsigned int i, unsigned int j) {
    int pos = rowPtr[i];
    bool found = 0;
    for (int n = rowPtr[i]; n < rowPtr[i] + (rowPtr[i + 1] - rowPtr[i]); n++) {
        if (colInd[n] == j) {
            pos = n;
            found = 1;
        }
    }
    if (found) return value[pos];
    return 0.0;
}

void CRSmatrix::changeValue(double x, unsigned int i, unsigned int j) {
    if (value.empty()) {
        value.push_back(x);
        colInd.push_back(j);
        rowPtr.push_back(0);
        rowPtr.push_back(1);
    } else {
        value.push_back(x);
        colInd.push_back(j);
        if (rowPtr.size() <= i + 1)
            rowPtr.push_back(value.size());
        else
            rowPtr[rowPtr.size() - 1] = rowPtr[rowPtr.size() - 1] + 1;
    }
}

std::vector<double> CRSmatrix::productAx(std::vector<double> x) {
    std::vector<double> product(x.size());
    for (int i = 0; i < rowNum; i++) {
        for (int j = 0; j < colNum; j++) {
            product[i] += x[j] * retrieveElement(i, j);
        }
    }
    return product;
}

void CRSmatrix::deleteValue(unsigned int i, unsigned int j) {
    for (int n = rowPtr[i]; n < rowPtr[i] + (rowPtr[i + 1] - rowPtr[i]); n++) {
        if (colInd[n] == j) {
            value.erase(value.begin() + n);
            colInd.erase(colInd.begin() + n);
            for (unsigned int m = i + 1; m < rowPtr.size(); m++) rowPtr[m]--;
        }
    }
}

void CRSmatrix::printA() {
    for (int valI = 0; valI < rowNum; valI++) {
        for (int valJ = 0; valJ < colNum; valJ++) {
            std::cout << " " << retrieveElement(valI, valJ);
        }
        std::cout << std::endl;
    }
}

double vectorNorm(std::vector<double> x) {
    double sum = 0.0;
    for (unsigned int i = 0; i < x.size(); i++) {
        sum += pow(x[i], 2);
    }
    return sqrt(sum);
}

std::vector<double> Jacobi(CRSmatrix A, std::vector<double> b) {
    //first find the inverse of the diagonal matrix
    CRSmatrix Di = A;
    for (int i = 0; i < Di.getRowSize(); i++) {
        for (int j = 0; j < Di.getColSize(); j++) {
            //remove all elements other than the diagonal
            if (i != j) {
                Di.deleteValue(i, j);
            } else {
                //take the inverse of the diagonal
                for (int n = Di.rowPtr[i]; n < Di.rowPtr[i] + (Di.rowPtr[i + 1] - Di.rowPtr[i]); n++) {
                    if (Di.colInd[n] == j) {
                        Di.value[n] = 1 / Di.value[n];
                    }
                }
            }
        }
    }

    //find L+U
    CRSmatrix LU = A;
    for (int i = 0; i < LU.getRowSize(); i++) {
        for (int j = 0; j < LU.getColSize(); j++) {
            if (i == j) {
                LU.deleteValue(i, j);
            } else {
                //take the inverse of the diagonal
                for (int n = LU.rowPtr[i]; n < LU.rowPtr[i] + (LU.rowPtr[i + 1] - LU.rowPtr[i]); n++) {
                    if (LU.colInd[n] == j) {
                        LU.value[n] = -1.0 * LU.value[n];
                    }
                }
            }
        }
    }

    //take D^-1 * (L+U), store in DLU
    CRSmatrix DLU = LU;
    for (int i = 0; i < DLU.getRowSize(); i++) {
        for (int j = 0; j < DLU.getColSize(); j++) {
            for (int n = DLU.rowPtr[i]; n < DLU.rowPtr[i] + (DLU.rowPtr[i + 1] - DLU.rowPtr[i]); n++) {
                if (DLU.colInd[n] == j) {
                    DLU.value[n] = Di.value[i] * LU.retrieveElement(i, j);
                }
            }
        }
    }

    //find intercept D^-1 * b
    std::vector<double> cep(b.size());
    cep = Di.productAx(b);

    //initial guess is the same as the intercept
    std::vector<double> x = cep;

    //find vetor norm ||b-Ax|| of initial guess
    std::vector<double> normTemp(x.size());
    normTemp = A.productAx(x);
    for (unsigned int i = 0; i < b.size(); i++) {
        normTemp[i] = -1 * normTemp[i];
        normTemp[i] = b[i] + normTemp[i];
    }
    int count = 0;

    //start Jacobi iteration
    std::vector<double> xTemp(b.size());
    std::vector<double> prod(b.size());
    double threshold = 10E-7;
    while (vectorNorm(normTemp) > threshold) {
        count++;

        //x(k-1)
        for (unsigned int i = 0; i < b.size(); i++) {
            xTemp[i] = x[i];
        }

        //find product of D^-1(L+U)*x
        prod = DLU.productAx(xTemp);

        //add D^-1(L+U)*x to D^-1 * b, store in x
        for (unsigned int i = 0; i < x.size(); i++) {
            x[i] = prod[i] + cep[i];
        }

        //find new vetor norm ||x(k-1)-x||
        for (unsigned int i = 0; i < b.size(); i++) {
            normTemp[i] = -1.0 * x[i];
            normTemp[i] = xTemp[i] + normTemp[i];
        }
    }
    return x;
}
