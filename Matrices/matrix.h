#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

using namespace std;

const bool DEBUG = false;

//Quick way to store a coefficient and its ratio to the scale value
struct zpair{
    double ratio = 0;
    int index = 0;
    zpair(double r, int i)
    {
        ratio = r;
        index = i;
    }
};


int getMaxIndex(vector<zpair> v);


template<class T>
void printV(vector<T> v)
{
    for(unsigned u = 0; u < v.size(); u++)
        cout << v.at(u) << " ";
    cout << endl;
}

//Class that represents a Matrix
class Matrix
{


public:
    /******************
     ** CONSTRUCTORS **
     ******************/
    Matrix();
    Matrix(int n, int m);
    Matrix(int n, int m, vector<double> meat);
    Matrix(int n, int m, vector<double> &meat, vector<double> &b);
    Matrix(Matrix &mat);
    Matrix(Matrix &mat, vector<double> &b);

    /***************
     ** ACCESSORS **
     ***************/
    void print();
    Matrix solve();

    Matrix operator=(Matrix &other);

private:

    int n = 0;                  //Number of rows
    int m = 0;                  //Number of columns
    int size = 0;               //m*n. Size of 1D array
    double* content = nullptr;  //Dynamically allocated 1D array that stores all
                                //values in the matrix
    vector<double> bValues;     //Answer vector

    /***************
     ** ACCESSORS **
     ***************/
    vector<double> getScale();
    double find_max(int start, int end);

    /***************
     ** MUTATORS  **
     ***************/
    void addRow(int cons, int change, double scalar);
    double getScalar(int cons, int change, int column);




};


#endif // MATRIX_H
