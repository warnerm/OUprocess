// C++ program to find Deteminant of a matrix. Adapted from http://www.geeksforgeeks.org/determinant-of-a-matrix/
#include "MHfuns.hpp"
// Function to get cofactor of mat[p][q] in temp[][]. n is current
// dimension of mat[][]
void getCofactor(double mat[nNode][nNode], double temp[nNode][nNode], int p, int q, int n)
{
    int i = 0, j = 0;
    
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = mat[row][col];
                
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
 n is current dimension of mat[][]. */
double determinantOfMatrix(double mat[nNode][nNode], int n)
{
    double D = 0; // Initialize result
    
    //  Base case : if matrix contains single element
    if (n == 1)
        return mat[0][0];
    
    double temp[nNode][nNode]; // To store cofactors
    
    int sign = 1;  // To store sign multiplier
    
    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of mat[0][f]
        getCofactor(mat, temp, 0, f, n);
        D += sign * mat[0][f] * determinantOfMatrix(temp, n - 1);
        
        // terms are to be added with alternate sign
        sign = -sign;
    }
    
    return D;
}


