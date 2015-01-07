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
      RefMat getViewLin(RefMat J, size_t i) const;
      const ConstRefMat getViewLin(const ConstRefMat J, size_t i) const;
      RefMat getViewNonLin(RefMat J, size_t i) const;
      const ConstRefMat getViewNonLin(const ConstRefMat J, size_t i) const;

      const Index& totalDimLin() const;
      const Index& totalDimNonLin() const;

    private:
      std::vector<Index> startLin_;
      std::vector<Index> dimLin_;
      Index totalDimLin_;
      std::vector<Index> startNonLin_;
      std::vector<Index> dimNonLin_;
      Index totalDimNonLin_;
  };
}


#endif // _PGS_CONSTRAINTMANAGER_H_
