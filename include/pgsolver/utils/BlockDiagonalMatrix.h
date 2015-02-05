#ifndef EIGEN_BLOCKDIAGONALMATRIX_H
#define EIGEN_BLOCKDIAGONALMATRIX_H

//assumption: square matrices, square blocks

namespace Eigen {

  namespace internal {

    template<typename Derived>
    class BlockDiagonalMatrixBase : public EigenBase<Derived>
    {
    public:
      enum {
        Flags = internal::traits<Derived>::Flags,
        CoeffReadCost = internal::traits<Derived>::CoeffReadCost,
        RowsAtCompileTime = internal::traits<Derived>::RowsAtCompileTime,
        ColsAtCompileTime = internal::traits<Derived>::ColsAtCompileTime,
        MaxRowsAtCompileTime = internal::traits<Derived>::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = internal::traits<Derived>::MaxColsAtCompileTime,
        MaxBlockSizeAtCompileTime = internal::traits<Derived>::MaxBlockSizeAtCompileTime,
        Options = internal::traits<Derived>::Options
      };
      typedef typename internal::traits<Derived>::Scalar Scalar;
      typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> DenseMatrixType;
      typedef typename DenseMatrixType::Index Index;
      typedef typename internal::traits<Derived>::CoefficientsType CoefficientsType;
      typedef typename internal::traits<Derived>::IndexVector IndexVector;
      typedef typename internal::traits<Derived>::BlockSizes BlockSizes;
      typedef EigenBase<Derived> Base;

    public:

      using Base::derived;
      using Base::rows;
      using Base::cols;

      /** \returns the number of blocks */
      inline Index blocks() const { return derived().blocks(); }

      /** \returns the \a i -th block (const version) */
      inline const Block<const CoefficientsType> block(Index i) const { return derived().block(i); }

      /** \returns the \a i -th block (non-const version) */
      inline Block<CoefficientsType> block(Index i) { return derived().block(i); }

      /** \return the \a i -th block as a fixed size block (const version)
        * \warning N must be the effective size of the block, i.e. equal to blockSize(i). */
      template<int N> inline const Block<const CoefficientsType, N, N> block(Index i) const { return derived().block<N>(i); }

      /** \return the \a i -th block as a fixed size block (non-const version)
        * \warning N must be the effective size of the block, i.e. equal to blockSize(i). */
      template<int N> inline Block<CoefficientsType, N, N> block(Index i) { return derived().block<N>(i); }

      /** \returns an expression of the underlying coefficient matrix */
      inline const CoefficientsType& coeffs() const { return derived().coeffs(); }

      /** \returns an expression of the underlying coefficient matrix */
      inline CoefficientsType& coeffs() { return derived().coeffs(); }

      /** \returns the vector of block sizes */
      inline const BlockSizes& blockSizes() const { return derived().blockSizes(); }

      /** \returns the vector of block start index */
      inline const IndexVector& blockStarts() const { return derived().starts(); }

      /** \returns the size of one block */
      inline Index blockSize(Index i) const 
      { 
        eigen_assert(0 <= i && i < blocks());
        return blockSizes()[i]; 
      }

      /** \returns the starting index of a block*/
      inline Index blockStart(Index i) const
      {
        eigen_assert(0 <= i && i < blocks());
        return blockStarts()[i];
      }

      template <typename Dest> inline void evalTo(Dest& dst) const
      {
        dst.resize(rows(), cols());
        dst.setZero();
        for (Index i = 0; i < blocks(); ++i)
        {
          Index start = blockStart(i);
          Index size = blockSize(i);
          dst.block(start, start, size, size) = block(i);
        }
      }

      DenseMatrixType toDenseMatrix() const
      {
        DenseMatrixType res(rows(), cols());
        evalTo(res);
        return res;
      }

    };



    //forward declaration
    template<typename, int, int, int, int, int> class BlockDiagonalMatrix;

    template<typename _Scalar, int _Size, int _BlocksAtCompileTime, int _MaxBlockSizeAtCompileTime, int _Options, int _MaxSizeAtCompileTime>
    struct traits<BlockDiagonalMatrix<_Scalar, _Size, _BlocksAtCompileTime, _MaxBlockSizeAtCompileTime, _Options, _MaxSizeAtCompileTime> >
    {
      typedef _Scalar Scalar;
      typedef Dense StorageKind;
      typedef DenseIndex Index;

      enum {
        Flags = LvalueBit,
        CoeffReadCost = NumTraits<Scalar>::ReadCost,
        RowsAtCompileTime = _Size,
        ColsAtCompileTime = _Size,
        MaxRowsAtCompileTime = _MaxSizeAtCompileTime,
        MaxColsAtCompileTime = _MaxSizeAtCompileTime,
        BlocksAtCompileTime = _BlocksAtCompileTime,
        MaxBlockSizeAtCompileTime = _MaxBlockSizeAtCompileTime,
        Options = _Options
      };
      typedef Matrix<Scalar, MaxBlockSizeAtCompileTime, ColsAtCompileTime, Options&RowMajor ? RowMajor : ColMajor, MaxBlockSizeAtCompileTime, MaxColsAtCompileTime> CoefficientsType;
      typedef Matrix<Index, BlocksAtCompileTime, 1, 0, (BlocksAtCompileTime != Dynamic) ? BlocksAtCompileTime : MaxRowsAtCompileTime, 1> IndexVector;
      typedef IndexVector BlockSizes;
    };

    template<typename _Scalar, int _Size=Dynamic, int _BlocksAtCompileTime=Dynamic, int _MaxBlockSizeAtCompileTime = Dynamic, int _Options = ColMajor, int _MaxSizeAtCompileTime = _Size>
    class BlockDiagonalMatrix : public BlockDiagonalMatrixBase<BlockDiagonalMatrix<_Scalar, _Size, _BlocksAtCompileTime, _MaxBlockSizeAtCompileTime, _Options, _MaxSizeAtCompileTime> >
    {
    public:
      typedef typename internal::traits<BlockDiagonalMatrix>::Scalar Scalar;
      typedef typename internal::traits<BlockDiagonalMatrix>::Index Index;
      typedef typename internal::traits<BlockDiagonalMatrix>::CoefficientsType CoefficientsType;
      typedef typename internal::traits<BlockDiagonalMatrix>::BlockSizes BlockSizes;
      typedef typename internal::traits<BlockDiagonalMatrix>::IndexVector IndexVector;

      enum {
        MaxBlockSizeAtCompileTime = internal::traits<BlockDiagonalMatrix>::MaxBlockSizeAtCompileTime
      };

      inline BlockDiagonalMatrix(const BlockSizes& blockSizes)
        : m_size(blockSizes.sum())
        , m_coeffs(((MaxBlockSizeAtCompileTime != Dynamic) ? 
              static_cast<Index>(MaxBlockSizeAtCompileTime) : m_size), m_size)
        , m_cumul(-1), m_nblocks(blockSizes.size())
        , m_blocksizes(blockSizes)
        , m_starts(blockSizes.size())
      {
        assert(_Size == Dynamic || m_size == _Size);
        eigen_assert(blockSizes.size() > 0);
        m_starts[0] = 0;
        for (Index i = 0; i < blockSizes.rows() - 1; ++i)   //[TODO]: we could implement a cumsum with compile-time optimization
          m_starts[i + 1] = m_starts[i] + blockSizes[i];
      }

      inline BlockDiagonalMatrix(Index size)
        : m_size(size)
        , m_coeffs(((MaxBlockSizeAtCompileTime != Dynamic) ? 
              static_cast<Index>(MaxBlockSizeAtCompileTime) : size), size)
        , m_cumul(0)
        , m_nblocks(0)
        , m_blocksizes(size)
        , m_starts(size)
      {
      }

      inline bool completed() const { return m_cumul == -1; }

      inline Block<CoefficientsType> push_block(Index blockSize)
      {
        eigen_assert(!completed() && "Block matrix is already complete");
        eigen_assert(blockSize <= MaxBlockSizeAtCompileTime || MaxBlockSizeAtCompileTime == Dynamic);
        eigen_assert(m_cumul + blockSize <= m_size);

        m_blocksizes[m_nblocks] = blockSize;
        m_starts[m_nblocks] = m_cumul;
        m_cumul += blockSize;
        ++m_nblocks;

        if (m_cumul == m_size)
          m_cumul = -1;         //complete

        return block_impl(m_nblocks - 1);
      }
      template<typename Derived>
      inline void push_block(const MatrixBase<Derived>& matrix)
      {
        assert(matrix.rows() == matrix.cols());
        Index blockSize = matrix.rows();
        push_block(blockSize) = matrix;
      }

      inline Index rows() const { return m_size; }
      inline Index cols() const { return m_size; }

      inline const BlockSizes& blockSizes() const { eigen_assert(completed()); return m_blocksizes; }
      inline const IndexVector& starts() const { eigen_assert(completed()); return m_starts; }
      inline Index blocks() const { eigen_assert(completed()); return m_nblocks; }

      inline const CoefficientsType& coeffs() const { return m_coeffs; }
      inline CoefficientsType& coeffs() { return m_coeffs; }

      inline const Block<const CoefficientsType> block(Index i) const 
      { 
        eigen_assert(completed());
        return m_coeffs.block(0, m_starts[i], m_blocksizes[i], m_blocksizes[i]); 
      }
      
      inline Block<CoefficientsType> block(Index i) 
      { 
        eigen_assert(completed());
        return block_impl(i);
      }

      template<int N> inline const Block<const CoefficientsType, N, N> block(Index i) const
      {
        eigen_assert(completed());
        eigen_assert(N == m_blocksizes[i]);
        return m_coeffs.template block<N, N>(0, m_starts[i]);
      }

      template<int N> inline Block<CoefficientsType, N, N> block(Index i)
      {
        eigen_assert(completed());
        eigen_assert(N == m_blocksizes[i]);
        return m_coeffs.template block<N, N>(0, m_starts[i]);
      }

      void resize(Index newSize)
      {
        m_coeffs.resize(((MaxBlockSizeAtCompileTime != Dynamic) ? 
              static_cast<Index>(MaxBlockSizeAtCompileTime) : m_size), m_size);
        m_cumul = 0;
        m_nblocks = 0;
        m_blocksizes.resize(newSize);
        m_starts.resize(newSize);
        m_size = newSize;
      }

    protected:
      inline Block<CoefficientsType> block_impl(Index i)
      {
        return m_coeffs.block(0, m_starts[i], m_blocksizes[i], m_blocksizes[i]);
      }


      Index m_size;
      CoefficientsType m_coeffs;
      Index m_cumul;
      Index m_nblocks;
      BlockSizes m_blocksizes;
      IndexVector m_starts;
    };

    template<typename Derived> class BlockDiagonalMatrixWrapper
      : public BlockDiagonalMatrixBase<BlockDiagonalMatrixWrapper<Derived> >
    {

    };
  }
}

#endif
