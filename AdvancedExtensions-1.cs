using System;
using System.Collections.Specialized;
using System.Globalization;
using System.Security.Cryptography.X509Certificates;
using Core;

namespace ProjectC {
    public static class AdvancedExtensions {
        /// <summary>
        /// This function creates the square submatrix given a square matrix as
        /// well as row and column indices to remove from it.
        /// </summary>
        ///
        /// <remarks>
        /// See page 246-247 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">An N-by-N matrix.</param>
        /// <param name="i">The index of the row to remove.</param>
        /// <param name="j">The index of the column to remove.</param>
        ///
        /// <returns>The resulting (N - 1)-by-(N - 1) submatrix.</returns>
        public static Matrix SquareSubMatrix(this Matrix a, int i, int j) {
            // A new matrix with a row and a column less is created
            Matrix newMatrix = new Matrix(a.M_Rows - 1, a.N_Cols - 1);

            // I have created two new variables to keep track of the rows and columns in the 
            // submatrix. This is to be able to skip a row/column in the original matrix, but keep
            // the row-/column-number the same in the submatrix 
            int newRow = 0;
            for (int x = 0; x < a.M_Rows; x++) {
                int newCol = 0;
                for (int y = 0; y < a.N_Cols; y++) {
                    // If x is the index of the row to remove x is incremented, so the row is 
                    // skipped 
                    if (x == i) {
                        x++;
                    }
                    // If the column is not the one that should be removed, and x is not the last
                    // row, the index is 
                    if (y != j && x != a.M_Rows) {
                        newMatrix[newRow, newCol] = a[x, y];
                        newCol++;
                    }
                }
                newRow++;
            }
            return newMatrix;
        }

        /// <summary>
        /// This function computes the determinant of a given square matrix.
        /// </summary>
        ///
        /// <remarks>
        /// See page 247 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <remarks>
        /// Hint: Use SquareSubMatrix.
        /// </remarks>
        ///
        /// <param name="a">An N-by-N matrix.</param>
        ///
        /// <returns>The determinant of the matrix</returns>
        public static double Determinant(this Matrix a) {
            double retval = 0;
            // Using recursion, if the matrix given contains just one index, that index is the 
            // determinant and thus is returned 
            if (a.M_Rows == 1) {
                retval = a[0, 0];
            } else {
                // Otherwise the determinant is calculated using definition 5.1 and 5.2 in the book 
                // page 247. The method SquareSubMatrix is then used to find the determinant of the 
                // submatrix (det(Aij) in definition 5.1. 
                for (int i = 1; i <= a.N_Cols; i++) {
                    // Goes through the top row in each (sub)matrix until the matrix is 1 x 1. 
                    retval += a[0, i - 1] * Math.Pow((-1), 1 + i) *
                              a.SquareSubMatrix(0, i - 1).Determinant
                                  ();
                }
            }
            return retval;
        }

        /// <summary>
        /// This function computes the Gram-Schmidt process on a given matrix.
        /// </summary>
        ///
        /// <remarks>
        /// See page 229 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">
        /// An M-by-N matrix. All columns are implicitly assumed linear
        /// independent.
        /// </param>
        ///
        /// <returns>
        /// A tuple (Q,R) where Q is a M-by-N orthonormal matrix and R is an
        /// N-by-N upper triangular matrix.
        /// </returns>
        public static Tuple<Matrix, Matrix> GramSchmidt(this Matrix a) {
            // We suggest this tolerance number for floating point comparisons.
            var tol = 1e-8;
            // The following is heavily inspired by the pseudocode algorithm found on page 229 of 
            // the book. 
            // Firstly, two matrices are created, one (Q) the size of a, the other (R) a square 
            // matrix with the same number of rows and columns as a has columns. 
            Matrix Q = new Matrix(a.M_Rows, a.N_Cols);
            Matrix R = new Matrix(a.N_Cols, a.N_Cols);

            // Q is set to be equal a. 
            for (int j = 0; j < a.N_Cols; j++) {
                for (int i = 0; i < a.M_Rows; i++) {
                    Q[i, j] = a[i, j];
                }
                // The dot product of the column q_i and the column u_i is found and saved in a 
                // variable called QU. 
                for (int i = 0; i < j; i++) {
                    double QU = 0;
                    for (int x = 0; x < a.M_Rows; x++) {
                        QU += Q[x, i] * a[x, j];
                    }
                    // r_ij is set to be QU, the dot product. 
                    R[i, j] = QU;
                    for (int y = 0; y < a.M_Rows; y++) {
                        Q[y, j] -= R[i, j] * Q[y, i];
                    }
                }
                // This bool and the loop following it is only used to decrease runtime, because if
                // the column q_j is a zero-column, we would waste our time calculating its  
                // magnitude. 
                bool jZero = true;
                for (int i = 0; i < a.M_Rows; i++) {
                    if (Math.Abs(Q[i, j]) > tol) {
                        jZero = false;
                        break;
                    }
                }
                if (jZero) {
                    break;
                } else {
                    // If the column is not zero, we calculate its magnitude and set r_jj to this 
                    // value. 
                    double norm = 0;
                    for (int i = 0; i < a.M_Rows; i++) {
                        norm += Math.Pow(Q[i, j], 2);
                    }

                    R[j, j] = Math.Sqrt(norm);
                    // Lastly q_j is set to q_j/r_jj
                    for (int i = 0; i < a.M_Rows; i++) {
                        Q[i, j] = Q[i, j] / R[j, j];
                    }
                }
            }
            return new Tuple<Matrix, Matrix>(Q, R);
        }
    }
}


