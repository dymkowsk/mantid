#ifndef MANTID_TESTV3D__
#define MANTID_TESTV3D__

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <ostream>
#include <vector>

#include "Matrix.h" 

using namespace Mantid;
using namespace Geometry;

class testMatrix: public CxxTest::TestSuite
{
private:


public:

  void makeMatrix(Matrix<double>& A) const
    {
      A.setMem(3,3);
      A[0][0]=1.0;
      A[1][0]=3.0;
      A[0][1]=4.0;
      A[0][2]=6.0;
      A[2][0]=5.0;
      A[1][1]=3.0;
      A[2][1]=1.0;
      A[1][2]=6.0;
      A[2][2]=-7.0;
      return;
    }

  /*!
    Test that a matrix can be inverted
  */
  void testInvert()
    {
      Matrix<double> A(3,3);
      
      A[0][0]=1.0;
      A[1][0]=3.0;
      A[0][1]=4.0;
      A[0][2]=6.0;
      A[2][0]=5.0;
      A[1][1]=3.0;
      A[2][1]=1.0;
      A[1][2]=6.0;
      A[2][2]=-7.0;

      Matrix<double> Iv=A;
      double Det=Iv.Invert();
      Matrix<double> Ident(3,3);
      Ident.identiyMatrix();
      
      TS_ASSERT(fabs(Det-105.0)<1e-5);
      TS_ASSERT(Ident==A);
    }

  /*!
    Check that we can swap rows and columns
  */
  void testSwapRows()
    {
      Matrix<double> A(3,3);
      makeMatrix(A);
      Matrix<double> B(A);
      A.swapRows(1,2);
      A.swapCols(1,2);
      TS_ASSERT_EQUALS(A[0][0],B[0][0]);
      TS_ASSERT_EQUALS(A[2][2],B[1][1]);
      // Plus all the others..
    }

  void testEigenvectors()
    {
      Matrix<double> Eval;
      Matrix<double> Diag;
      Matrix<double> A(3,3);  // NOTE: A must be symmetric
      A[0][0]=1.0;
      A[1][0]=A[0][1]=4.0;
      A[0][2]=A[2][0]=5.0;
      A[1][1]=3.0;
      A[2][1]=A[1][2]=6.0;
      A[2][2]=-7.0;
      TS_ASSERT(A.Diagonalise(Eval,Diag));

      Matrix<double> MA=A*Eval;
      Matrix<double> MV=Eval*Diag;
      Eval.sortEigen(Diag);
      TS_ASSERT(Diag[0][0]<Diag[1][1]);
      TS_ASSERT(Diag[1][1]<Diag[2][2]);
      TS_ASSERT(MA==MV);

      std::vector<double> X(3);
      X[0]=Eval[0][1];
      X[1]=Eval[1][1];
      X[2]=Eval[2][1];

      std::vector<double> out=A*X;
      transform(X.begin(),X.end(),X.begin(),std::bind2nd(std::multiplies<double>(),Diag[1][1]));
      TS_ASSERT_EQUALS(X[0],out[0]);
      TS_ASSERT_EQUALS(X[1],out[1]);
      TS_ASSERT_EQUALS(X[2],out[2]);

    }
  
  /*!
    Tests the diagonalisation  on a symmetric 2x2 matrix
  */
  void testDiagonalise()
    {
      Matrix<double> Eval;
      Matrix<double> Diag;
      Matrix<double> A(2,2);   // symmetric only
      A[0][0]=1.0;
      A[1][0]=3.0;
      A[0][1]=3.0;
      A[1][1]=4.0;
      TS_ASSERT(A.Diagonalise(Eval,Diag));  // returns 1 or 2

      Matrix<double> EvalT(Eval);
      EvalT.Transpose();
      Eval*=Diag;
      Eval*=EvalT;
      TS_ASSERT(Eval==A);
    }

}

	
};

#endif
