using System;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Security.Policy;
using Core;

namespace ProjectB
{
    public static class GaussExtensions
    {
        /// <summary>
        /// This function computes the elementary row replacement operation on
        /// the given matrix.
        /// </summary>
        ///
        /// <remarks>
        /// Note that we add the row (as in the lectures) instead of subtracting
        /// the row (as in the textbook).
        /// </remarks>
        ///
        /// <param name="a">
        /// An N-by-M matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">
        /// The index of the row to replace.
        /// </param>
        /// <param name="m">
        /// The multiplum of row j to add to row i.
        /// </param>
        /// <param name="j">
        /// The index of the row to replace with.
        /// </param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after having performed the elementary
        /// row operation.
        /// </returns>
        public static Matrix ElementaryRowReplacement(
            this Matrix a, int i, double m, int j) {
            // Counts up in the i'th row to the number of columns-1 and adds another rows 
            // corresponding column multiplied with a double. 
            for (int r = 0; r < a.N_Cols; r++) {
                a[i, r] =  m * a[j, r] + a[i, r]; 
            }
            return a;
        }

        /// <summary>
        /// This function computes the elementary row interchange operation on
        /// the given matrix.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">
        /// The index of the first row of the rows to interchange.
        /// </param>
        /// <param name="j">
        /// The index of the second row of the rows to interchange.
        /// </param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after having performed the elementary
        /// row operation.
        /// </returns>
        public static Matrix ElementaryRowInterchange(
            this Matrix a, int i, int j) {
            for (int r = 0; r < a.N_Cols; r++) {
                // temp is a variable that is used, so the first rows indexes are not lost 
                var temp = a[j, r];
                a[j, r] = a[i, r];
                a[i, r] = temp;

            }

            return a; 
        }

        /// <summary>
        /// This function computes the elementary row scaling operation on the
        /// given matrix.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">The index of the row to scale.</param>
        /// <param name="c">The value to scale the row by.</param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after having performed the elementary
        /// row operation.
        /// </returns>
        public static Matrix ElementaryRowScaling(
            this Matrix a, int i, double c) {
            for (int r = 0; r < a.N_Cols; r++) {
                a[i, r] = a[i, r] * c;
            }

            return a;
        }


        /// <summary>
        /// This function executes the forward reduction algorithm provided in
        /// the assignment text to achieve row Echelon form of a given
        /// augmented matrix.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M augmented matrix.
        /// </param>
        ///
        /// <returns>
        /// An N-by-M matrix that is the row Echelon form.
        /// </returns>
        public static Matrix ForwardReduction(this Matrix a)
        {
            // We suggest this tolerance number for floating point comparisons.
            var tol = 1e-8;
            var pivotRow = 0; 

            // This is against the suggestion to iterate over rows in the outer-loop, but it is 
            // more important to find the row in which there is a non-zero entry than the column. 
            for (int y = 0; y < a.N_Cols; y++) {
                for (int x = pivotRow; x < a.M_Rows; x++) {
                    // tol is used to compare instead of 0, as suggested above. 
                    if (a[x, y] > tol || -tol > a[x, y]) {
                        // Switches the non-zero-index row with the lowest row that has not been 
                        // interchanged yet. 
                        a.ElementaryRowInterchange(x, pivotRow);
                        // Increments pivotRow, because the first row that should be checked is now
                        // one row lower
                        pivotRow++;
                        // Perform ElementaryRowReplacement on all the indexes in the rows below the
                        // one just found, so all entries below the pivotindex are 0. 
                        for (int i = pivotRow; i < a.M_Rows; i++) {
                            a.ElementaryRowReplacement(i, -a[i, y]/a[pivotRow-1,y], pivotRow-1);
                        }
                        break;
                    }
                }
            }
            return a; 
        }

        /// <summary>
        /// This function executes the backward reduction algorithm provided in
        /// the assignment text given an augmented matrix in row Echelon form.
        /// </summary>
        ///
        /// <param name="a">
        /// An N-by-M augmented matrix in row Echelon form.
        /// </param>
        ///
        /// <returns>
        /// The resulting N-by-M matrix after executing the algorithm.
        /// </returns>
        public static Matrix BackwardReduction(this Matrix a) {
            // We suggest this tolerance number for floating point comparisons.
            var tol = 1e-8;
            var pivotColumn = a.N_Cols;
            var pivotRow = a.M_Rows;

            // This loop follows the suggestion to iterate over rows in the outer loop, but does 
            // it backwards, so it starts from the bottom row
            for (int i = pivotRow-1; i > -1; i--) {
                // The column loop starts from the column furthest to the left 
                for (int j = 0; j < pivotColumn; j++) {
                    var tempIndex = a[i, j];
                    if (tempIndex > tol || tempIndex < -tol) {
                        a.ElementaryRowScaling(i, 1 / tempIndex);
                        // After the row has been scaled, every entry directly above the pivot entry
                        // is RowReplaced, so it contains a 0. 
                        for (int x = i-1; x > -1; x--) {
                            a.ElementaryRowReplacement(x, -a[x,j], i);
                        }
                        break;
                    }
                }
            }
            return a;
        }
        

        /// <summary>
        /// This function performs Gauss elimination of a linear system
        /// given in matrix form by a coefficient matrix and a right hand side
        /// vector. It is assumed that the corresponding linear system is
        /// consistent and has exactly one solution.
        /// </summary>
        ///
        /// <remarks>
        /// Hint: Combine ForwardReduction and BackwardReduction.
        /// </remarks>
        ///
        /// <param name="a">An N-by-M matrix.</param>
        /// <param name="b">An N-size vector.</param>
        ///
        /// <returns>The M-sized vector x such that a * x = b.</returns>
        public static Vector GaussElimination(this Matrix a, Vector b) {

            // Creates a new matrix consisting of the matrix a on the right and the vector b on the 
            // left
            Matrix newMatrix = a.AugmentRight(b);

            newMatrix.ForwardReduction();
            
            newMatrix.BackwardReduction();
            
            // After the newmatrix has a reduced row echelon form a new vector is created 
            Vector x = new Vector(a.N_Cols);

            // The new vector consists of the entries in the column furthest to the left
            for (int l = 0; l < a.M_Rows; l++) {
                x[l] = newMatrix[l, newMatrix.N_Cols-1];
            }
            return x;


        }
    }

    public static class BasicExtensions
    {
        /// <summary>
        /// This function creates an augmented matrix given a matrix 'a' and a
        /// right-hand side vector 'v'.
        /// </summary>
        ///
        /// <remarks>
        /// See page 12 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">An M-by-N matrix.</param>
        /// <param name="v">An M-size vector.</param>
        ///
        /// <returns>The M-by-(N + 1) augmented matrix [a | v].</returns>
        public static Matrix AugmentRight(this Matrix a, Vector v)
        {
            var mRows = a.M_Rows;
            var nCols = a.N_Cols;

            var retval = new double[mRows, nCols + 1]; // 0-initialized

            for (var i = 0; i < mRows; i++)
            {
                for (var j = 0; j < nCols; j++)
                {
                    retval[i, j] = a[i, j];
                }
                retval[i, nCols] = v[i];
            }

            return new Matrix(retval);
        }
    }
}
