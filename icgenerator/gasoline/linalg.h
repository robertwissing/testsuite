
/** linalg.h
 * 
 *  PKDGRAV Source Code
 * 
 *  Author: Kenneth W. Flynn
 *          flynnk@astro.umd.edu
 *  Mods:   Derek C. Richardson
 *          dcr@astro.umd.edu
 *          [3-D hardwired for speed]
 *
 *  Linear algebra support routines for PKDGRAV.  Used initially for 
 *   aggregates, but of more general use.  In case of question, these routines
 *   are based on the Numerical Recipes methods, but are original work in
 *   terms of source code.  The primary reason for this is the handling of
 *   vectors and matrices in the NRiC functions, particularly the non-zero
 *   indexing.  Rather than copy the matrices, pass them in, and then have to
 *   copy them back to the more traditional form, we've simply implemented the
 *   routines ourselves.
 *
 *  Output parameters need to be allocated, but not initialized.
 */
 
#ifndef __LINALG_H
#define __LINALG_H
 
#include "floattype.h"
 
/* maximum number of Jacobi sweeps for routine jacobi() */
#define MAX_JACOBI_SWEEPS 50

typedef FLOAT Scalar;
typedef Scalar Vector[3];
typedef Vector Matrix[3]; /* can't declare const Matrix with this defn (why?) */

/** Copy a vector.
 *
 *  Parameters (in):
 *   u - The source vector
 *
 *  Parameters (out):
 *   v - The destination vector
 */
void vectorCopy(const Vector u,Vector v);

/** Multiply the specified vector times a scalar and store the result in
 *   destination.  
 *
 *  Parameters (in):
 *   u - The vector to multiply
 *   scalar - Scalar to multiply by
 *
 *  Parameters (out):
 *   v - Destination vector for the result.  May be same as u if desired
 */
void vectorScale(const Vector u,FLOAT scalar,Vector v);

/** Add two vectors together and store the result in the specified location.
 *
 *  Parameters (in):
 *   v1 - First vector to add
 *   v2 - Second vector to add
 *
 *  Parameters (out):
 *   v - Destination for sum, may be same as either v1 or v2 if desired
 */
void vectorAdd(const Vector v1,const Vector v2,Vector v);

/** Subtract one vector from a second and store the result in the specified location.
 *
 *  Parameters (in):
 *   v1 - Vector to subtract from
 *   v2 - Vector to subtract
 *
 *  Parameters (out):
 *   v - Destination for difference, may be same as either v1 or v2 if desired
 */
void vectorSub(const Vector v1,const Vector v2,Vector v);

/** Returns dot product of two vectors.
 *
 *  Parameters (in):
 *   v1 - First vector
 *   v2 - Second vector
 *
 *  Returns:
 *   The dot product
 */
FLOAT vectorDot(const Vector v1,const Vector v2);

/*DEBUG add comments*/
FLOAT vectorMagSq(const Vector v);
FLOAT vectorMag(const Vector v);
void vectorNorm(Vector v);

/** Computes cross product of two vectors.
 *
 *  Parameters (in):
 *   v1 - First vector
 *   v2 - Second vector
 *
 *  Parameters (out):
 *   v - Destination for the cross product
 */
void vectorCross(const Vector v1,const Vector v2,Vector v);

/*DEBUG add comment*/
void vectorSet(Vector v,double x,double y,double z);

/** Sets the specified vector's components to zero.
 *
 *  Parameters (out):
 *   v - Zero vector
 */
void vectorZero(Vector v);

/*DEBUG add comment*/
void vectorGetBasis(Vector a,Vector b,Vector c);

/** Copy the specified matrix into the second matrix.
 *
 *  Parameters (in):
 *   a - The source matrix
 * 
 *  Parameters (out):
 *   b - Destination matrix
 */
void matrixCopy(Matrix a,Matrix b);

/** Multiply the given matrix times the given vector and store the result
 *   in the supplied destination.
 *
 *  Parameters (in):
 *   m - The matrix to multiply by
 *   u - The vector being multiplied
 *
 *  Parameters (out):
 *   v - Place to store the results, which must be different from u
 *        in order to generate correct results
 */
void matrixTransform(Matrix m,const Vector u,Vector v);

/** Take the transpose of a matrix.
 *
 *  Parameters (in):
 *   a - The matrix in question
 *
 *  Parameters (out):
 *   b - Its transpose
 */
void matrixTranspose(Matrix a,Matrix b);

/** Swap two rows in a matrix. 
 *
 *  Parameters (in/out):
 *   m - The source matrix, which is modified to contain resulting matrix
 *
 *  Parameters (in):
 *   row1 - First row
 *   row2 - Second row
 */
void matrixSwapRows(Matrix m,int row1,int row2);

/** Generate identity matrix.
 *
 *  Parameters (out):
 *   m - Identity matrix
 */
void matrixIdentity(Matrix m);

/** Generate diagonal matrix where each element, m[i][i], equals
 *   the corresponding value in the specified vector, v[i].
 *
 *  Parameters (in):
 *   v - Vector of values for the diagonal of the matrix
 *
 *  Parameters (out):
 *   m - Diagonal matrix fill with the specified values
 */
void matrixDiagonal(const Vector v,Matrix m);

/** Sum the off diagonal elements of the specified matrix.
 *
 *  Parameters:
 *   m - Matrix in question
 *
 *  Returns:
 *   The sum of any elements not on the diagonal of the matrix.
 */ 
FLOAT matrixSumOffDiagElem(Matrix m);

/** Sum the absolute value of the off diagonal elements of the specified
 *   matrix.
 *
 *  Parameters:
 *   m - Matrix in question
 *
 *  Returns:
 *   The sum of the absolute value of any elements not on the diagonal of 
 *    the matrix.
 */
FLOAT MatrixSumAbsOffDiagElem(Matrix m);

/** Multiply the specified matrix times a scalar and store the result in
 *   destination.  
 *
 *  Parameters (in):
 *   a - The matrix to multiply
 *   s - Scalar to multiply by
 *
 *  Parameters (out):
 *   b - Destination matrix for the result.  May be same as a if desired
 */
void matrixScale(Matrix a,Scalar s,Matrix b);

/** Multiply two matrices and store the result in the specified matrix. 
 *
 *  Parameters (in):
 *   a - First matrix
 *   b - Second matrix
 * 
 *  Parameters (out):
 *   c - a x b
 */
void matrixMultiply(Matrix a,Matrix b,Matrix c);

/** Take the inverse of a matrix, using Gauss-Jordan elimination with partial
 *   pivoting.
 *
 *  Parameters (in):
 *   mat_in - The matrix in question
 *
 *  Parameters (out):
 *   mat_out - Inverse
 */
void matrixInverse(Matrix mat_in,Matrix mat_out);

/** Compute the eigenvalues and eigenvectors of a real symmetric matrix
 *   using the Jacobi method.
 * 
 *  Parameters (in):
 *   m - A real symmetric matrix
 *
 *  Parameters (out):
 *   eig_vals - The 3 eigenvalues of the system
 *   eig_vecs - Matrix whose columns are the eigenvectors of the system
 */
void jacobi(Matrix m,Vector eig_vals,Matrix eig_vecs); 

#endif

void InvertSymMatrix3x3( 
    double m11, double m12, double m13, 
    double m22, double m23, double m33, 
    double im11, double im12, double im13,
    double im22, double im23, double im33 ) ;

/* Note: Using det!=0 as estimator for good behaviour
  If determinant non-zero not sufficient, could use 
  condition number to detect pathological cases 
  Condition number can be estimated as CN = |m|.|minv|
  with Frobenius norm ||=sqrt(sum m_ij^2) */
#define InvertSymMatrix3x3( m11, m12, m13, m22, m23, m33, im11, im12, im13, im22, im23, im33 ) { \
    double det = m11* (m22* m33- m23* m23)- \
             m12* (m12* m33- m23* m13)+ \
             m13* (m12* m23- m22* m13); \
    double idet; \
    if(det > 0.0) {idet = 1/det;} else{idet = 0.0;} \
    im11= (m22* m33- m23* m23)* idet; \
    im12= (m13* m23- m12* m33)* idet; \
    im13= (m12* m23- m13* m22)* idet; \
    im22= (m11* m33- m13* m13)* idet; \
    im23= (m12* m13- m11* m23)* idet; \
    im33= (m11* m22- m12* m12)* idet; }

void InvertMatrix3x3( 
    double m11, double m12, double m13, 
    double m21, double m22, double m23, 
    double m31, double m32, double m33, 
    double im11, double im12, double im13, 
    double im21, double im22, double im23, 
    double im31, double im32, double im33 );

#define InvertMatrix3x3( m11, m12, m13, m21, m22, m23, m31, m32, m33, im11, im12, im13, im21, im22, im23, im31, im32, im33 ) { \
    double idet, det = m11* (m22* m33- m32* m23)- \
             m12* (m21* m33- m23* m31)+ \
             m13* (m21* m32- m22* m31); \
    assert(det != 0); \
    idet = 1/det; \
    im11= (m22* m33- m32* m23)* idet; \
    im12= (m13* m32- m12* m33)* idet; \
    im13= (m12* m23- m13* m22)* idet; \
    im21= (m23* m31- m21* m33)* idet; \
    im22= (m11* m33- m13* m31)* idet; \
    im23= (m21* m13- m11* m23)* idet; \
    im31= (m21* m32- m31* m22)* idet; \
    im32= (m31* m12- m11* m32)* idet; \
    im33= (m11* m22- m21* m12)* idet; }

void FrobeniusNorm(
    double m11, double m12, double m13,
    double m22, double m23, double m33) ;

#define FrobeniusNorm(m11, m12, m13, m22, m23, m33) (sqrt((m11)*(m11) + 2*(m12)*(m12) + 2*(m13)*(m13) + (m22)*(m22) + 2*(m23)*(m23) + (m33)*(m33)))
