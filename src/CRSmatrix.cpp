/*
 *  CRSmatrix.cpp
 *  CRS matrix functions and Jacobi method
 *
 */

#include "CRSmatrix.h"

using namespace std;

CRSmatrix::CRSmatrix(vector<double> v, vector<int> r, vector<int> c)
{
    value = v;
    rowPtr = r;
    colInd = c;

    //find total number of rows
    rowNum = rowPtr.size()-1;

    //find total number of columns
    colNum = 0;
    for (int i = 0; i <= colInd.size(); i++) {
    if (colInd[i]>colNum) colNum = colInd[i];
    }
    colNum++;
}

CRSmatrix::CRSmatrix(int r, int c)
{
    rowNum = r;
    colNum = c;
}

CRSmatrix::CRSmatrix(string mtxFilePath) {
    //Read in files
    std::ifstream mtxFile(mtxFilePath);

    // Check if the file is open
    if (!mtxFile) {
        cout << "Unable to open file" << endl;
        exit(1);
    }

    // Ignore headers and comments:
    while (mtxFile.peek() == '%') mtxFile.ignore(2048, '\n');

    // Read defining parameters:
    int rowRank, columnRank, nNZ;
    mtxFile >> rowRank >> columnRank >> nNZ;

    while(mtxFile.good()){
        for (int i = 0; i < nNZ; i++)
        {
            int mtxRow, mtxCol;
            double mtxValue;
            mtxFile >> mtxRow >> mtxCol >> mtxValue;
            rowPtr.push_back(mtxRow-1); // row and col are 1 based in mtx files
            colInd.push_back(mtxCol-1);
            value.push_back(mtxValue);
            cout << mtxRow-1 << " " << mtxCol-1 << " " << mtxValue << endl;
        }
    }

    //find total number of rows
    rowNum = rowPtr.size()-1;

    //find total number of columns
    colNum = 0;
    for (int i = 0; i < colInd.size(); i++) {
        if (colInd[i] > colNum) colNum = colInd[i];
    }
    colNum++;

    cout << rowNum << " " << colNum << endl;
}

CRSmatrix::CRSmatrix(string valueAddress, string rowPtrAddress, string colIndAddress) {
    //Read in files
    ifstream rowFile, valueFile, colFile;

    rowFile.open(rowPtrAddress);
    valueFile.open(valueAddress);
    colFile.open(colIndAddress);

    // Check if the file is open
    if (!rowFile || !valueFile || !colFile) {
        cout << "Unable to open file" << endl;
        exit(1);
    }

    // Declare a variable to load the contents from the file
    string line = "";

    // Loading the value vector
    string myValue;
    while(valueFile.good()){
        getline(valueFile,myValue,'\n');
        //convert string value to double value
        double myDoubleValue = atof(myValue.c_str());
        value.push_back(myDoubleValue);
    }

    // Loading the row pointer vector
    while (rowFile >> line) {
        int rowNum = stoi(line);
        rowPtr.push_back(rowNum-1);
    }

    // Loading the column indice vector
    while (colFile >> line) {
        int colIndex = stoi(line);
        colInd.push_back(colIndex-1);
    }

    //find total number of rows
    rowNum = rowPtr.size()-1;

    //find total number of columns
    colNum = 0;
    for (int i = 0; i < colInd.size(); i++) {
        if (colInd[i]>colNum) colNum = colInd[i];
    }
    colNum++;
    cout << colNum << endl;
}

double CRSmatrix::retrieveElement (int i, int j) {
    int pos = rowPtr[i];
    bool found = 0;
    for (int n = rowPtr[i]; n<rowPtr[i]+ (rowPtr[i+1] - rowPtr[i]); n++){
        if (colInd[n] == j){
            pos = n;
            found = 1;
        }
    }
    if (found) return value[pos];
    return 0.0;
}

void CRSmatrix::changeValue (double x, int i, int j) {
    bool found = 0;
    if (value.empty()){
        value.push_back(x);
        colInd.push_back(j);
        rowPtr.push_back(0);
        rowPtr.push_back(1);
    } else {
        value.push_back(x);
        colInd.push_back(j);
        if (rowPtr.size() <= i+1) rowPtr.push_back(value.size());
        else rowPtr[rowPtr.size()-1] = rowPtr[rowPtr.size()-1]+1;
    }
}

vector<double> CRSmatrix:: productAx(vector<double> x){
    vector<double> product(x.size());
    for (int i = 0; i<rowNum; i++){
        for (int j = 0; j<colNum; j++){
            product[i] += x[j]*retrieveElement(i, j);
        }
    }
    return product;
}

void CRSmatrix::deleteValue(int i, int j){
    for (int n = rowPtr[i]; n<rowPtr[i]+ (rowPtr[i+1] - rowPtr[i]); n++){
        if (colInd[n] == j){
            value.erase(value.begin()+n);
            colInd.erase(colInd.begin()+n);
            for (int m = i+1; m<rowPtr.size(); m++) rowPtr[m]--;
        }
    }
}

void CRSmatrix::printA(){
    for (int valI = 0; valI<rowNum; valI++){
        for (int valJ = 0; valJ<colNum; valJ++){
            cout<<" "<<retrieveElement(valI, valJ);
        }
        cout<<endl;
    }
}

double vectorNorm(vector<double> x) {
    double sum = 0.0;
    for (int i = 0; i<x.size(); i++){
        sum += pow(x[i], 2);
    }
    return sqrt(sum);
}

vector<double> Jacobi(CRSmatrix A, vector<double> b) {
    //first find the inverse of the diagonal matrix
    CRSmatrix Di = A;
    for (int i = 0; i<Di.rowNum; i++){
        for (int j = 0; j<Di.colNum; j++){
            //remove all elements other than the diagonal
            if (i != j){
                Di.deleteValue(i, j);
            }else{
                //take the inverse of the diagonal
                int pos;
                for (int n = Di.rowPtr[i]; n<Di.rowPtr[i]+(Di.rowPtr[i+1] - Di.rowPtr[i]); n++){
                    if (Di.colInd[n] == j){
                        pos = n;
                        Di.value[pos] = 1/Di.value[pos];
                    }
                }
            }
        }
    }

    //find L+U
    CRSmatrix LU = A;
    for (int i = 0; i<LU.rowNum; i++) {
        for (int j = 0; j < LU.colNum; j++) {
            if (i == j) {
                LU.deleteValue(i, j);
            }else {
                //take the inverse of the diagonal
                int pos;
                for (int n = LU.rowPtr[i]; n < LU.rowPtr[i] + (LU.rowPtr[i + 1] - LU.rowPtr[i]); n++) {
                    if (LU.colInd[n] == j) {
                        pos = n;
                        LU.value[pos] = -1.0 * LU.value[pos];
                    }
                }
            }
        }
    }

    //take D^-1 * (L+U), store in DLU
    CRSmatrix DLU = LU;
    for (int i = 0; i<DLU.rowNum; i++){
        for (int j = 0; j<DLU.colNum; j++){
            for (int n = DLU.rowPtr[i]; n<DLU.rowPtr[i]+(DLU.rowPtr[i+1] - DLU.rowPtr[i]); n++){
                if (DLU.colInd[n] == j){
                    DLU.value[n] = Di.value[i]*LU.retrieveElement(i, j);
                }
            }
        }
    }

    //find intercept D^-1 * b
    vector<double> cep(b.size());
    cep = Di.productAx(b);

    //initial guess is the same as the intercept
    vector<double> x = cep;

    //find vetor norm ||b-Ax|| of initial guess
    vector<double> normTemp(x.size());
    normTemp = A.productAx(x);
    for (int i = 0; i < b.size(); i++) {
        normTemp[i] = -1*normTemp[i];
        normTemp[i] = b[i] + normTemp[i];
    }
    int count=0;

    //start Jacobi iteration
    vector<double> xTemp(b.size());
    vector<double> prod (b.size());
    while (vectorNorm(normTemp) > 10E-7){
        count++;

        //x(k-1)
        for (int i = 0; i < b.size(); i++) {
            xTemp[i] = x[i];
        }

        //find product of D^-1(L+U)*x
        prod = DLU.productAx(xTemp);

        //add D^-1(L+U)*x to D^-1 * b, store in x
        for (int i = 0; i<x.size(); i++){
            x[i] = prod[i] + cep[i];
        }

        //find new vetor norm ||x(k-1)-x||
        for (int i = 0; i < b.size(); i++) {
            normTemp[i] = -1.0*x[i];
            normTemp[i] = xTemp[i] + normTemp[i];
        }
    }
    return x;
}
