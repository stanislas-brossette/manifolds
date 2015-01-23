#ifndef _PGS_CONSTRAINTMANAGER_H_
#define _PGS_CONSTRAINTMANAGER_H_

#include <vector>

#include <Eigen/Core>
#include <pgsolver/manifolds/defs.h>
#include <pgsolver/solver/Problem.h>

namespace pgs
{
  class ConstraintManager
  {
    public:
      ConstraintManager();
      void init(Problem& problem);
      /// \brief returns the part of J corresponding to the Linear Part of
      /// constraint J with the assumption that J contains only linear parts of
      /// constraints
      // TODO: find a cleaner method to differentiate the const and nonConst
      // versions of the method
      RefMat getViewLin(RefMat J, size_t i) const;
      const ConstRefMat getConstViewLin(const ConstRefMat J, size_t i) const;
      /// \brief returns the part of J corresponding to the NonLinear Part of
      /// constraint J with the assumption that J contains only Nonlinear parts of
      /// constraints
      // TODO: find a cleaner method to differentiate the const and nonConst
      // versions of the method
      RefMat getViewNonLin(RefMat J, size_t i) const;
      const ConstRefMat getConstViewNonLin(const ConstRefMat J, size_t i) const;

      const Index& totalDimLin() const;
      const Index& totalDimNonLin() const;
      const Index& totalDim() const;

    private:
      std::vector<Index> startLin_;
      std::vector<Index> dimLin_;
      Index totalDimLin_;
      std::vector<Index> startNonLin_;
      std::vector<Index> dimNonLin_;
      Index totalDimNonLin_;
      Index totalDim_;
  };
}


#endif // _PGS_CONSTRAINTMANAGER_H_
