#include <iostream>
#include "matrix.h"
using namespace std;

void printV(vector<double> v);
int main()
{
    vector<double> content = { 3, -13, 9,   3,
                              -6,   4, 1, -18,
                               6,  -2, 2,   4,
                              12,  -8, 6,  10};
    Matrix m(4 ,4 ,content);

    vector<double> b = {-19, -34, 16, 26};
    m.solve(b);
    return 0;
}

void printV(vector<double> v)
{
    for(unsigned u = 0; u < v.size(); u++)
        cout << v.at(u) << " ";
    cout << endl;
}
