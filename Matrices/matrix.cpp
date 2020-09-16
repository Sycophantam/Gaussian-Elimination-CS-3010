#include "matrix.h"
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <cmath>

Matrix::Matrix()
{
    content = new double(0);
    size = n = m = 1;
}
Matrix::Matrix(int n,   //Number of rows in the matrix
               int m)   //Number of columns in the matrix
{
    this->n = n;
    this->m = m;
    //Creates a 2D array to represent a n x m matrix
    //This will be annoying to get a certain value, but it's cheaper with memory
    content = new double[n * m];
    for (int i = 0; i < n * m; i++)
        *(content + i) = 0;
    size = n*m;
}
Matrix::Matrix(int n,               //Number of rows in the matrix
               int m,               //Number of columns in the matrix
               vector<double> meat) //Values to be put in the matrix
{
    this->n = n;
    this->m = m;
    size = n*m;

    //Checking if the vector's content can evenly fit in the matrix
    if(meat.size() != size)
    {
        cout << "Error: Vector of size " << meat.size() << " cannot fit in " <<
                n << " by " << m << " matrix" << endl;
        size = this->n = this->m = 0;
        return;
    }
    content = new double[size];
    std::copy(meat.begin(), meat.end(), content);


}

Matrix::Matrix(Matrix &mat)
{
    n = mat.n;
    m = mat.m;
    size = n*m;


}

void Matrix::print(vector<double> &b)
{
    //We want to make sure there is content before we print the matrix
    assert(size != 0);
    for(int i = 0; i < n; i++)
    {
        cout << "[";
        for(int j = 0; j < m; j++)
        {
            cout << setw(7) << left << setprecision(2) << fixed << *(content + (i*m + j));
        }
        cout << setw(7) << right << setprecision(2) << fixed << b.at(i);
        cout << "]" << endl;
    }

}

Matrix Matrix::solve(vector<double> b)  //Vector of b values
{
    //Variable that is stored as the last value in the exclusion vector

    print(b);
    //Scale vector
    vector<double> scale = getScale();
    if(DEBUG)
    {
        cout << "Scale: ";
        for(unsigned i = 0; i < scale.size(); i++)
        {
            cout << scale.at(i) << " ";
        }
        cout << endl;
    }

    //Vector containing the excluded indices. Starts off with bad due to the way
    //c++ finds values
    vector<int> exc;

    //Vector containing the coefficients
    vector<double> coef;

    //Vector containing the ratios and the indeces the ratio comes from
    vector<zpair> pivot;
    //Traversing the matrix horizontally
    for(int i = 0; i < m; i++)
    {
        //Getting the coefficients of each column
        for(int j = 0; j < n; j++)
        {
            if(!std::count(exc.begin(), exc.end(), j))
            {
                double add = *(content + i + j*m);
                coef.push_back(add);

                //Storing the ratio and the index where it came from
                pivot.push_back(zpair(abs(add)/scale.at(j), j));
            }
            if(DEBUG)
            {
                cout << "Exclusion vector:";
                printV(exc);
                cout << "Pivot indices: ";
                for(int afa = 0; afa < pivot.size(); afa++)
                {
                    cout << pivot.at(afa).index << " ";
                }
                cout << endl;
                cout << "Pivot ratios: ";
                for(int afa = 0; afa < pivot.size(); afa++)
                {
                    cout << pivot.at(afa).ratio << " ";
                }
                cout << endl;

            }
        }
        cout << "Pivot ratios: ";
        for(int afa = 0; afa < pivot.size(); afa++)
        {
            cout << pivot.at(afa).ratio << " ";
        }
        cout << endl;
        int pivotRow = getMaxIndex(pivot);
        cout << "Pivot row is " << pivotRow << endl;
        cout << endl;


        //Excluding the largest row from being counted again. We add it to the
        //front so that bad is constantly at the back
        exc.insert(exc.begin(), pivotRow);
        //printV(exc);
        int scale = 0;
        for(int k = 0; k < n; k++)
        {
            if(DEBUG)
            {
                cout << "In for loop" << endl;
                cout << "exc.size() is " << exc.size() << endl;
                cout << "k is " << k << endl;
                cout << "exc is ";
                printV(exc);
                cout << "count is ";
                cout << std::count(exc.begin(), exc.end(), k) << endl;
            }
            //If the index is not excluded, we reduce the row
            if(!std::count(exc.begin(), exc.end(), k))
            {
                if(DEBUG)
                {
                    cout << "If statement passed" << endl;
                }
                double scalar = getScalar(pivotRow, k, i);
                if(DEBUG)
                {
                    cout << "Scalar is " << scalar << endl;
                    cout << "k is " << k << endl;
                    cout << "pivotRow is " << pivotRow << endl;
                }
//                print(b);
//                cout << endl;
                addRow(pivotRow, k, scalar, b);
                print(b);
                cout << endl;
            }
        }
        pivot.clear();
        coef.clear();
    }
}

double Matrix::getScalar(int cons, int change, int column)
{
    if(DEBUG)
    {
        cout << "Getting the scalar" << endl;
        cout << "The numerator should be " << *(content + change*m + column);
        cout << endl;
        cout << "The denominator should be " << *(content + cons*m + column);
        cout << endl;
    }
    return -1 * ((*(content + change*m + column))/
                 (*(content + cons*m + column)));
}
void Matrix::addRow(int cons,           //Row that is adding to change
                    int change,         //Row that is being added to
                    double scalar,      //Scale factor of the row
                    vector<double> &b)  //Answer vector
{
    for(int i = 0; i < m; i++)
    {
        //Change*m + i gives the index that is being added to
        //Cons*m + i gives the index that is adding to change
        *(content + change*m + i) += scalar * (*(content + cons*m + i));
    }
    b.at(change) += scalar * b.at(cons);
}

//Finds the index of the maximum value of the given vector
int getMaxIndex(vector<zpair> v)
{
    double max = -1000000;
    unsigned index = 0;
    for(unsigned i = 0; i < v.size(); i++)
    {
        if(v.at(i).ratio > max)
        {
            max = v.at(i).ratio;
            index = i;
        }
    }
    return static_cast<int>(v.at(index).index);
}
vector<double> Matrix::getScale()
{
    vector<double> scale;
    for(int i = 0; i < n; i++)
    {
        scale.push_back(find_max(i*m, i*m + m));
    }
    return scale;
}

double Matrix::find_max(int start, int end)
{
    double max = -100000;
    for(int i = start; i < end; i++)
    {
        if(abs(*(content + i)) > max)
            max = *(content + i);
    }
    return abs(max);
}

Matrix Matrix::operator=(Matrix &other)
{
    Matrix m;
    m.n = other.n;
    m.m = other.m;
    m.size = m.n*m.m;

    //Checking if the vector's content can evenly fit in the matrix
    m.content = new double[size];
    std::copy(other.getVec().begin(), other.getVec().end(), m.content);
    return m;
}

vector<double> Matrix::getVec()
{
    return vector<double> (content, content + size);
}
