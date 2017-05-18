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

// Function to get adjoint of A[nNode][nNode] in adj[nNode][nNode].
//From http://www.geeksforgeeks.org/adjoint-inverse-matrix/
void adjoint(double A[nNode][nNode])
{
    // temp is used to store cofactors of A[][]
    int sign = 1;
    double temp[nNode][nNode];
    
    for (int i=0; i<nNode; i++)
    {
        for (int j=0; j<nNode; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, nNode);
            
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
            
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinantOfMatrix(temp, nNode-1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
void inverse(double A[nNode][nNode])
{
    // Find determinant of A[][]
    double det = determinantOfMatrix(A, nNode);
    
    // Find adjoint
    adjoint(A);
    
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<nNode; i++)
        for (int j=0; j<nNode; j++)
            inv[i][j] = adj[i][j]/(det);
}


