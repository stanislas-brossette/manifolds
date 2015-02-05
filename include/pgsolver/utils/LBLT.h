#ifndef EIGEN_LBLT_H
#define EIGEN_LBLT_H

namespace Eigen {

  namespace internal {
    template<typename MatrixType, int UpLo> struct LBLT_Traits;
  }

  template<typename _MatrixType, int _UpLo> class LBLT
  {
  public:
    typedef _MatrixType MatrixType;
    enum {
      RowsAtCompileTime = MatrixType::RowsAtCompileTime,
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      Options = MatrixType::Options & ~RowMajorBit, // these are the options for the TmpMatrixType, we need a ColMajor matrix here!
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
      UpLo = _UpLo
    };
    typedef typename MatrixType::Scalar Scalar;
    typedef typename NumTraits<typename MatrixType::Scalar>::Real RealScalar;
    typedef typename MatrixType::Index Index;
    typedef Matrix<Scalar, RowsAtCompileTime, 2, Options, MaxRowsAtCompileTime, 2> TmpMatrixType;

    typedef Transpositions<RowsAtCompileTime, MaxRowsAtCompileTime, Index> TranspositionType;

    typedef internal::LBLT_Traits<MatrixType, UpLo> Traits;
    typedef typename Traits::MatrixB MatrixB;

    /** \brief Default Constructor.
    *
    * The default constructor is useful in cases in which the user intends to
    * perform decompositions via LBLT::compute(const MatrixType&).
    */
    LBLT()
      : m_matrix(),
      m_transpositions(),
      m_isInitialized(false)
    {}

    /** \brief Default Constructor with memory preallocation
    *
    * Like the default constructor but with preallocation of the internal data
    * according to the specified problem \a size.
    * \sa LBLT()
    */
    LBLT(Index size)
      : m_matrix(size, size),
      m_blocks(size),
      m_transpositions(size),
      m_temporary(size,2),
      m_isInitialized(false)
    {}

    /** \brief Constructor with decomposition
    *
    * This calculates the decomposition for the input \a matrix.
    * \sa LBLT(Index size)
    */
    LBLT(const MatrixType& matrix)
      : m_matrix(matrix.rows(), matrix.cols()),
      m_blocks(matrix.rows()),
      m_transpositions(matrix.rows()),
      m_temporary(matrix.rows(),2),
      m_isInitialized(false)
    {
      compute(matrix);
    }

    /** \returns a view of the upper triangular matrix U */
    inline typename Traits::MatrixU matrixU() const
    {
      eigen_assert(m_isInitialized && "LBLT is not initialized.");
      return Traits::getU(m_matrix);
    }

    /** \returns a view of the lower triangular matrix L */
    inline typename Traits::MatrixL matrixL() const
    {
      eigen_assert(m_isInitialized && "LBLT is not initialized.");
      return Traits::getL(m_matrix);
    }

    /** \returns the permutation matrix P as a transposition sequence.*/
    inline const TranspositionType& transpositionsP() const
    {
      eigen_assert(m_isInitialized && "LBLT is not initialized.");
      return m_transpositions;
    }

    /** \returns the coefficients of the block diagonal matrix B */
    inline const typename Traits::MatrixB& matrixB() const
    {
      eigen_assert(m_isInitialized && "LBLT is not initialized.");
      return m_blocks;
    }

    /** \returns a solution x of \f$ A x = b \f$ using the current decomposition of A.
    *
    * This function also supports in-place solves using the syntax <tt>x = decompositionObject.solve(x)</tt> .
    *
    * \note_about_checking_solutions
    *
    * More precisely, this method solves \f$ A x = b \f$ using the decomposition \f$ A = P^T L B L^* P \f$
    * by solving the systems \f$ P^T y_1 = b \f$, \f$ L y_2 = y_1 \f$, \f$ B y_3 = y_2 \f$,
    * \f$ L^* y_4 = y_3 \f$ and \f$ P x = y_4 \f$ in succession. If the matrix \f$ A \f$ is singular, then
    * \f$ B \f$ will also be singular (all the other matrices are invertible). In that case, the
    * least-square solution of \f$ B y_3 = y_2 \f$ is computed. This does not mean that this function
    * computes the least-square solution of \f$ A x = b \f$ is \f$ A \f$ is singular.
    */
    template<typename Rhs>
    inline const internal::solve_retval<LBLT, Rhs>
      solve(const MatrixBase<Rhs>& b) const
    {
        eigen_assert(m_isInitialized && "LBLT is not initialized.");
        eigen_assert(m_matrix.rows() == b.rows()
          && "LBLT::solve(): invalid number of rows of the right hand side matrix b");
        return internal::solve_retval<LBLT, Rhs>(*this, b.derived());
      }

    template<typename Derived>
    bool solveInPlace(MatrixBase<Derived> &bAndX) const;

    LBLT& compute(const MatrixType& matrix);

    MatrixType reconstructedMatrix() const;

    inline Index rows() const { return m_matrix.rows(); }
    inline Index cols() const { return m_matrix.cols(); }

    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "LDLT is not initialized.");
      return Success;
    }

  protected:

    MatrixType m_matrix;
    MatrixB m_blocks;
    TranspositionType m_transpositions;
    TmpMatrixType m_temporary;
    //internal::SignMatrix m_sign;
    bool m_isInitialized;
  };

  namespace internal {

    template<int UpLo> struct lblt_inplace;

    template<> struct lblt_inplace<Lower>
    {
      template<typename MatrixType, typename BlockMatrixType, typename TranspositionType, typename Workspace>
      static bool unblocked(MatrixType& mat, BlockMatrixType& B, TranspositionType& transpositions, Workspace& temp)
      {
        using std::abs;
        using std::sqrt;
        typedef typename MatrixType::Scalar Scalar;
        typedef typename MatrixType::RealScalar RealScalar;
        typedef typename MatrixType::Index Index;
        eigen_assert(mat.rows() == mat.cols());
        const Index size = mat.rows();

        const double alpha = (1. + sqrt(17.)) / 8.;
        Index k = 0;
        while (k < size - 1)
        {
          Index s = 1;
          //find the largest off-diagonal element in column k
          Index r;
          RealScalar lambda = mat.col(k).tail(size - k - 1).cwiseAbs().maxCoeff(&r);
          r += k+1;

          if (lambda>0)
          {
            if (abs(mat(k, k)) >= alpha*lambda)
            {
              s = 1;
              transpositions.coeffRef(k) = k;
            }
            else
            {
              //find the largest off-diagonal element of in row r. 
              //since we read only the lowest part of mat, we need to cut this computation in two
              Index p, p1, p2;
              RealScalar sigma;
              RealScalar sigma1 = mat.row(r).segment(k, r - k).cwiseAbs().maxCoeff(&p1);
              RealScalar sigma2 = mat.col(r).tail(size - r - 1).cwiseAbs().maxCoeff(&p2);
              if (sigma1 > sigma2)
              {
                sigma = sigma1; p = p1 + k;
              }
              else
              {
                sigma = sigma2; p = p2 + r + 1;
              }

              if (sigma*abs(mat(k, k)) >= alpha*lambda)
              {
                s = 1;
                transpositions.coeffRef(k) = k;
              }
              else if (abs(mat(r, r)) >= alpha*sigma)
              {
                s = 1;
                transpositions.coeffRef(k) = r;
              }
              else
              {
                s = 2;
                if (p == k)
                {
                  if (r == k + 1)
                  {
                    transpositions.coeffRef(k) = k; transpositions.coeffRef(k+1) = k + 1;
                  }
                  else
                  {
                    transpositions.coeffRef(k) = k; transpositions.coeffRef(k+1) = r;
                  }
                }
                else if (p == k + 1)
                {
                  transpositions.coeffRef(k) = r; transpositions.coeffRef(k+1) = r;
                }
                else
                {
                  if (r == k + 1)
                  {
                    transpositions.coeffRef(k) = p; transpositions.coeffRef(k+1) = k + 1;
                  }
                  else
                  {
                    transpositions.coeffRef(k) = p; transpositions.coeffRef(k+1) = r;
                  }
                }
              }
            }

            //std::cout << "%before pivot " << k << std::endl;
            //std::cout << "M = " << (toMatlab)mat << std::endl << std::endl;
            //applying transposition while taking care of manipulating only the lowest part
            for (Index i = 0; i < s; ++i)
            {
              Index piv = transpositions.coeff(k + i);
              if (piv != k + i)
              {
                Index ki = k + i;
                Index t = size - piv - 1; // trailing size after the biggest element
                mat.col(ki).tail(t).swap(mat.col(piv).tail(t));
                mat.row(ki).head(ki).swap(mat.row(piv).head(ki));
                std::swap(mat.coeffRef(ki, ki), mat.coeffRef(piv,piv));
                for (Index j = ki + 1; j<piv; ++j)
                {
                  Scalar tmp = mat.coeffRef(j, ki);
                  mat.coeffRef(j, ki) = numext::conj(mat.coeffRef(piv, j));
                  mat.coeffRef(piv, j) = numext::conj(tmp);
                }
                if (NumTraits<Scalar>::IsComplex)
                  mat.coeffRef(piv, ki) = numext::conj(mat.coeff(piv, ki));
              }
            }

            //std::cout << "%before update " << k << std::endl;
            //std::cout << "M = " << (toMatlab)mat << std::endl << std::endl;
            Index t = size - k - s;
            SelfAdjointView<Block<MatrixType, Dynamic, Dynamic>, Lower> C = mat.bottomRightCorner(t, t).template selfadjointView<Lower>();
            if (s == 1)
            {
              Scalar Akk = mat.coeffRef(k, k);
              RealScalar invAkk = 1. / numext::real(Akk);
              C.rankUpdate(mat.col(k).tail(t), -invAkk);
              //C.rankUpdate(mat.block<Dynamic, 1>(k+1,k,t,1), -invAkk);
              //for (Index j = k + 1; j < size; ++j)
              //{
              //  //mat.col(j).tail(size - j) -= invAkk*mat.coeffRef(j, k)* mat.col(k).tail(size - j);
              //  mat.block<Dynamic,1>(j,j,size-j,1) -= invAkk*mat.coeffRef(j, k)* mat.block<Dynamic,1>(j,k,size-j,1);           
              //}
              mat.col(k).tail(t) *= invAkk;
              B.push_block(1).coeffRef(0,0) = Akk;
            }
            else //s == 2
            {
              Matrix<Scalar,2,2> D = mat.template block<2, 2>(k, k).template selfadjointView<Lower>();
              Matrix<Scalar, 2, 2> invD = D.inverse();
              C.rankUpdate(mat.col(k).tail(t), -invD(0,0));
              C.rankUpdate(mat.col(k+1).tail(t), -invD(1,1));
              C.rankUpdate(mat.col(k).tail(t), mat.col(k+1).tail(t), -invD(1,0));
              temp.bottomRows(t) = mat.template block<Dynamic, 2>(k + 2, k, t, 2) * invD;
              mat.template block<Dynamic, 2>(k + 2, k, t, 2) = temp.bottomRows(t);
              mat(k + 1, k) = Scalar(0);
              B.push_block(D);
            }
          }
          else
          {
            transpositions.coeffRef(k) = k;
            B.push_block(1).coeffRef(0, 0) = mat.coeffRef(k, k);
          }
          //std::cout << "%end of loop " << k  << std::endl;
          //std::cout << "M = " << (toMatlab)mat << std::endl << std::endl;
          k += s;
        }

        //if there remains a 1x1 bottom right corner to process 
        //(i.e. exit from the main loop after a step with s=1 or for matrix of size 1)
        if (k == size - 1)
        {
          transpositions.coeffRef(k) = k;
          B.push_block(1).coeffRef(0, 0) = mat.coeffRef(k, k);
        }
        return true;
      }
    };

    //template<> struct lblt_inplace<Upper>
    //{
    //  template<typename MatrixType, typename TranspositionType, typename Workspace>
    //  static EIGEN_STRONG_INLINE bool unblocked(MatrixType& mat, TranspositionType& transpositions, Workspace& temp, SignMatrix& sign)
    //  {
    //    Transpose<MatrixType> matt(mat);
    //    return lblt_inplace<Lower>::unblocked(matt, transpositions, temp, sign);
    //  }
    //};

    template<typename MatrixType> struct LBLT_Traits<MatrixType, Lower>
    {
      typedef const TriangularView<const MatrixType, UnitLower> MatrixL;
      typedef const TriangularView<const typename MatrixType::AdjointReturnType, UnitUpper> MatrixU;
      typedef BlockDiagonalMatrix<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime,
                  Dynamic, 2, MatrixType::Options, MatrixType::MaxRowsAtCompileTime> MatrixB;
      static inline MatrixL getL(const MatrixType& m) { return m; }
      static inline MatrixU getU(const MatrixType& m) { return m.adjoint(); }
    };

    template<typename MatrixType> struct LBLT_Traits<MatrixType, Upper>
    {
      typedef const TriangularView<const typename MatrixType::AdjointReturnType, UnitLower> MatrixL;
      typedef const TriangularView<const MatrixType, UnitUpper> MatrixU;
      typedef BlockDiagonalMatrix<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime,
                  Dynamic, 2, MatrixType::Options, MatrixType::MaxRowsAtCompileTime> MatrixB;
      static inline MatrixL getL(const MatrixType& m) { return m.adjoint(); }
      static inline MatrixU getU(const MatrixType& m) { return m; }
    };
  } // end namespace internal

  /** Compute / recompute the LBLT decomposition A = L D L^* = U^* D U of \a matrix
    */
  template<typename MatrixType, int _UpLo>
  LBLT<MatrixType, _UpLo>& LBLT<MatrixType, _UpLo>::compute(const MatrixType& a)
  {
    eigen_assert(a.rows() == a.cols());
    const Index size = a.rows();

    m_matrix = a;

    m_transpositions.resize(static_cast<int>(size));
    m_isInitialized = false;
    m_temporary.resize(size,2);
    m_blocks.resize(size);

    internal::lblt_inplace<UpLo>::unblocked(m_matrix, m_blocks, m_transpositions, m_temporary);

    m_isInitialized = true;
    return *this;
  }

  /** \returns the matrix represented by the decomposition,
  * i.e., it returns the product: P^T L B L^* P.
  * This function is provided for debug purpose. */
  template<typename MatrixType, int _UpLo>
  MatrixType LBLT<MatrixType, _UpLo>::reconstructedMatrix() const
  {
    eigen_assert(m_isInitialized && "LBLT is not initialized.");
    const Index size = m_matrix.rows();
    MatrixType res(size, size);

    // P
    res.setIdentity();
    res = transpositionsP() * res;
    // L^* P
    res = matrixU() * res;
    // B(L^*P)
    for (Index i = 0; i < m_blocks.blocks(); ++i)
    {
      res.middleRows(m_blocks.blockStart(i), m_blocks.blockSize(i)) = m_blocks.block(i)*res.middleRows(m_blocks.blockStart(i), m_blocks.blockSize(i));
    }
    // L(BL^*P)
    res = matrixL() * res;
    // P^T (LBL^*P)
    res = transpositionsP().transpose() * res;

    return res;
  }
}

#endif
