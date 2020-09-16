#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

using namespace std;

const bool DEBUG = false;

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
class Matrix
{


public:
    Matrix();
    Matrix(int n, int m);
    Matrix(int n, int m, vector<double> meat);
    Matrix(Matrix &mat);
    void print(vector<double> &b);
    Matrix solve(vector<double> b);
    void addRow(int cons, int change, double scalar, vector<double> &b);
    Matrix operator=(Matrix &other);
    vector<double> getVec();

private:

    int n = 0;
    int m = 0;
    int size = 0;
    double* content = nullptr;
    vector<double> getScale();
    double getScalar(int cons, int change, int column);
    double find_max(int start, int end);



};


#endif // MATRIX_H
