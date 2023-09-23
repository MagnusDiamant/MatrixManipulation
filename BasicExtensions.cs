using System;
using System.Security.Cryptography;
using Core;

namespace ProjectA
{
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

        /// <summary>
        /// This function computes the matrix-vector product of a matrix 'a' and
        /// a column vector 'v'.
        /// </summary>
        ///
        /// <remarks>
        /// See page 68 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">An M-by-N matrix.</param>
        /// <param name="v">An N-size vector.</param>
        ///
        /// <returns>The M-sized vector a * v.</returns>
        public static Vector Product(this Matrix a, Vector v) {
            var mRows = a.M_Rows;
            var nCols = a.N_Cols;
            // A check to see if the matrix and the vector are compatible ie. the number of rows 
            // in the matrix is equivalent to the number of values in the vector. 
            if (nCols == v.Size) {
                
                // 
                var retval = (new double[mRows]);
                // The value on the i'th row and j'th column in the matrix is multiplied with the 
                // j'th value in the vector. The product is added to retval. 
                for (var i = 0; i < mRows; i++) {
                    for (var j = 0; j < nCols; j++) {
                        retval[i] += a[i, j] * v[j];
                    }
                }

                return new Vector(retval);
            } 
            // If the matrix and vector are not compatible, an exception is thrown. 
            else {
                throw new ArgumentException("Matrix's number of columns and Vector's" +
                                            " number of rows must be the same");
            }
            
        }

        /// <summary>
        /// This function computes the matrix product of two given matrices 'a'
        /// and 'b'.
        /// </summary>
        ///
        /// <remarks>
        /// See page 58 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">An M-by-N matrix.</param>
        /// <param name="b">An N-by-P matrix.</param>
        ///
        /// <returns>The M-by-P matrix a * b.</returns>
        public static Matrix Product(this Matrix a, Matrix b) {
            var amRows = a.M_Rows;
            var bnCols = b.N_Cols;
            var aCols = a.N_Cols;
            // Checks if the matrices are compatible, ie. matrix a's number of columns is the same 
            // as matrix b's number of rows 
            if (aCols == b.M_Rows) {
                var retval = (new double[amRows, bnCols]);

                // Uses three for-loops to keep track. The i-loop keeps track of the rows in the 
                // product-matrix and the rows in the a-matrix. The x-loop keeps track of the 
                // columns in the product-matrix and the columns in the b-matrix. Lastly, the 
                // j-loop keeps track of the columns in the a-matrix and the rows in the b-matrix. 
                for (var i = 0; i < amRows; i++) {
                    for (var x = 0; x < bnCols; x++) {
                        for (var j = 0; j < aCols; j++) {
                            retval[i, x] += a[i, j] * b[j, x];
                        }
                    }
                }
                return new Matrix(retval);
            } 
            // If the two matrices are not compatible, and exception is thrown. 
            else {
                throw new ArgumentException("Matrix a's number of columns and matrix b's" +
                                            " number of rows must be the same");
            }
        }

        /// <summary>
        /// This function computes the transpose of a given matrix.
        /// </summary>
        ///
        /// <remarks>
        /// See page 69 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="a">An M-by-N matrix.</param>
        ///
        /// <returns>The N-by-M matrix a^T.</returns>
        public static Matrix Transpose(this Matrix a) {
            // The variables' names are changes to represent the rows and columns of the transposed
            // matrix
            var aCols = a.M_Rows;
            var aRows = a.N_Cols;

            var retval = (new double[aRows, aCols]);

            // This for-loop flips the matrix across the diagonal
            for (var i = 0; i < aRows; i++) {
                for (var j = 0; j < aCols; j++) {
                    retval[i, j] = a[j, i];
                }
            }
            return new Matrix(retval);
        }

        /// <summary>
        /// This function computes the Euclidean vector norm of a given vector.
        /// </summary>
        ///
        /// <remarks>
        /// See page 197 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="v">An N-dimensional vector.</param>
        ///
        /// <returns>The Euclidean norm of the vector.</returns>
        public static double VectorNorm(this Vector v)
        {
            var norm = new double();
            // This for-loop takes the power of each value in the vector and adds it to the 
            // variable norm
            for (var i = 0; i < v.Size; i++) {
                norm += Math.Pow(v[i], 2); 
            }
            // Returns the square-root of norm
            return Math.Sqrt(norm);
        }
    }
}
