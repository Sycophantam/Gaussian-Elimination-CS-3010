#include "matrix.h"
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <cmath>

/**********************************************************
 *
 * Constructor Matrix: Class Matrix
 *_________________________________________________________
 * This method constructs the Matrix class
 *_________________________________________________________
 * PRE-CONDITIONS
 *     none
 *
 * POST-CONDITIONS
 *     This function will construct a Matrix class
 ***********************************************************/
Matrix::Matrix()
{
    content = new double(0);
    size = n = m = 1;
}

/**********************************************************
 *
 * Constructor Matrix: Class Matrix
 *_________________________________________________________
 * This method constructs the Matrix class
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      n(int)-Number of rows in the matrix
 *      m(int)-Number of columns in the matrix
 *
 * POST-CONDITIONS
 *     This function will construct a Matrix class
 ***********************************************************/
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

/**********************************************************
 *
 * Constructor Matrix: Class Matrix
 *_________________________________________________________
 * This method constructs the Matrix class
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      n(int)-Number of rows in the matrix
 *      m(int)-Number of columns in the matrix
 *      meat(vector<double>)-Values in the matrix
 *
 * POST-CONDITIONS
 *     This function will construct a Matrix class
 ***********************************************************/
Matrix::Matrix(int n,               //Number of rows in the matrix
               int m,               //Number of columns in the matrix
               vector<double> meat) //Values to be put in the matrix
{
    this->n = n;
    this->m = m;
    size = n*m;

    //Checking if the vector's content can evenly fit in the matrix
    if(meat.size() != static_cast<unsigned>(size))
    {
        cout << "Error: Vector of size " << meat.size() << " cannot fit in " <<
                n << " by " << m << " matrix" << endl;
        size = this->n = this->m = 0;
        return;
    }
    content = new double[size];
    std::copy(meat.begin(), meat.end(), content);


}

/**********************************************************
 *
 * Constructor Matrix: Class Matrix
 *_________________________________________________________
 * This method constructs the Matrix class
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      n(int)-Number of rows in the matrix
 *      m(int)-Number of columns in the matrix
 *      meat(vector<double>)-Values in the matrix
 *      b(vector<double>)-Vector of answers
 *
 * POST-CONDITIONS
 *     This function will construct a Matrix class
 ***********************************************************/
Matrix::Matrix(int n,                   //Number of rows in the matrix
               int m,                   //Number of columns in the matrix
               vector<double> &meat,    //Values to be put in the matrix
               vector<double> &b)       //Vector of answers
{
    Matrix mat(n, m, meat);
    *this = mat;
    bValues = b;

}

/**********************************************************
 *
 * Constructor Matrix: Class Matrix
 *_________________________________________________________
 * This method constructs the Matrix class
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      mat(Matrix)-Matrix to copy
 *
 * POST-CONDITIONS
 *     This function will construct a Matrix class
 ***********************************************************/
Matrix::Matrix(Matrix &mat)     //Matrix to copy
{
    n = mat.n;
    m = mat.m;
    size = n*m;
    if(content != nullptr)
    {
        delete content;
        content = nullptr;
    }

    content = new double[size];

    //Copying the values from the original matrix into the copied
    for(int i = 0; i < mat.size; i++)
    {
        *(content + i) = *(mat.content + i);
    }
    bValues = mat.bValues;

}

/**********************************************************
 *
 * Constructor Matrix: Class Matrix
 *_________________________________________________________
 * This method constructs the Matrix class
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      mat(Matrix)-Matrix to copy
 *      b(vector<double>)-Answer vector
 *
 * POST-CONDITIONS
 *     This function will construct a Matrix class
 ***********************************************************/
Matrix::Matrix(Matrix &mat,         //Matrix to copy
               vector<double> &b)   //Vector as answers
{
    *this = mat;
    bValues = b;
}

/**********************************************************
 *
 * Method print(): Class Matrix
 *_________________________________________________________
 * This method prints the matrix and its values
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      None
 *
 * POST-CONDITIONS
 *     This function returns nothing
 ***********************************************************/
void Matrix::print()
{
    //We want to make sure there is content before we print the matrix
    assert(size != 0);
    for(int i = 0; i < n; i++)
    {
        cout << "[";
        for(int j = 0; j < m; j++)
        {
            cout << setw(7) << left
                 << setprecision(2) << fixed << *(content + (i*m + j));
        }
        if(!bValues.empty())
            cout << setw(7) << right
                 << setprecision(2) << fixed << bValues.at(i);
        cout << "]" << endl;
    }

}

/**********************************************************
 *
 * Method solve(): Class Matrix
 *_________________________________________________________
 * This method uses Gaussian elimination with scaled partial
 * pivoting to solve the system of matrix
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      None
 *
 * POST-CONDITIONS
 *     This function returns the modified matrix
 ***********************************************************/
Matrix Matrix::solve()  //Vector of b values
{
    //Variable that is stored as the last value in the exclusion vector

    print();
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
                for(unsigned afa = 0; afa < pivot.size(); afa++)
                {
                    cout << pivot.at(afa).index << " ";
                }
                cout << endl;
                cout << "Pivot ratios: ";
                for(unsigned afa = 0; afa < pivot.size(); afa++)
                {
                    cout << pivot.at(afa).ratio << " ";
                }
                cout << endl;

            }
        }
        cout << "Pivot ratios: ";
        for(unsigned u = 0; u < pivot.size(); u++)
        {
            cout << pivot.at(u).ratio << " ";
        }
        cout << endl;
        int pivotRow = getMaxIndex(pivot);
        cout << "Pivot row is " << pivotRow + 1 << endl;
        cout << endl;


        //Excluding the largest row from being counted again. We add it to the
        //front so that bad is constantly at the back
        exc.insert(exc.begin(), pivotRow);
        //printV(exc);
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
                addRow(pivotRow, k, scalar);
                print();
                cout << endl;
            }
        }
        pivot.clear();
        coef.clear();
    }
    cout << "Exclusion vector" << endl;
    printV(exc);
    return *this;
}

/**********************************************************
 *
 * Method getScalar(): Class Matrix
 *_________________________________________________________
 * This method finds the factor needed to reduce a row
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      cons(int)-The row that is not changed
 *      change(int)-The row that is modified
 *      column(int)-The column to reduce
 *
 * POST-CONDITIONS
 *     This function returns the scalar needed to reduce the first nonzero value
 *      in the row to zero
 ***********************************************************/
double Matrix::getScalar(int cons,      //Row that will not be modified
                         int change,    //Row that will be modified
                         int column)    //Column to base the scalar off of
{
    return -1 * ((*(content + change*m + column))/
                 (*(content + cons*m + column)));
}

/**********************************************************
 *
 * Method addRow(): Class Matrix
 *_________________________________________________________
 * This method adds row cons to row change
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      cons(int)-The row that is adding
 *      change(int)-The row that is changed
 *      scalar(double)-The factor needed to reduce one term to zero when adding
 *                      cons to change
 *
 * POST-CONDITIONS
 *     This function returns nothing
 ***********************************************************/
void Matrix::addRow(int cons,           //Row that is adding to change
                    int change,         //Row that is being added to
                    double scalar)      //Scale factor of the row
{
    for(int i = 0; i < m; i++)
    {
        //Change*m + i gives the index that is being added to
        //Cons*m + i gives the index that is adding to change
        *(content + change*m + i) += scalar * (*(content + cons*m + i));
    }
    if(!bValues.empty())
        bValues.at(change) += scalar * bValues.at(cons);
}

/**********************************************************
 *
 * Method getMaxIndex(): Class Matrix
 *_________________________________________________________
 * This method finds the index of the maximum ratio of coefficient to scale
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      v(zpair)-Vector of ratios and indices where they come from
 *
 * POST-CONDITIONS
 *     This function returns the index of the largest ratio
 ***********************************************************/
int getMaxIndex(vector<zpair> v)    //Vector containing the ratios and the
                                    //index where they came from
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

/**********************************************************
 *
 * Method getScale(): Class Matrix
 *_________________________________________________________
 * This method finds the scale vector
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      None
 *
 * POST-CONDITIONS
 *     This function returns the scale vector
 ***********************************************************/
vector<double> Matrix::getScale()
{
    vector<double> scale;
    for(int i = 0; i < n; i++)
    {
        scale.push_back(find_max(i*m, i*m + m));
    }
    return scale;
}

/**********************************************************
 *
 * Method find_max(): Class Matrix
 *_________________________________________________________
 * This method finds the maximum value in each row of the matrix
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      start(int)-Starting location to look for
 *      end(int)-Ending location of the list
 *
 * POST-CONDITIONS
 *     This function returns the maximum value between content[start]
 *      and content[end]
 ***********************************************************/
double Matrix::find_max(int start,  //Index of content to start looking
                        int end)    //Index of content to stop looking
{
    double max = -100000;
    for(int i = start; i < end; i++)
    {
        if(abs(*(content + i)) > max)
            max = *(content + i);
    }
    return abs(max);
}

/**********************************************************
 *
 * Operator =: Class Matrix
 *_________________________________________________________
 * Overloaded assignment operator
 *_________________________________________________________
 * PRE-CONDITIONS
 *     The following need pre-defined values:
 *      other(Matrix)-Matrix to copy
 *
 * POST-CONDITIONS
 *     This function returns a copy of other
 ***********************************************************/
Matrix Matrix::operator=(Matrix &other)
{
    n = other.n;
    m = other.m;
    size = n*m;
    bValues = other.bValues;

    if(content != nullptr)
    {
        delete content;
        content = nullptr;
    }

    content = new double[size];
    for(int i = 0; i < other.size; i++)
    {
        *(content + i) = *(other.content + i);
    }
    return *this;
}

