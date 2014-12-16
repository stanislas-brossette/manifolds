#ifndef _PGS_CONSTRAINTMANAGER_H_
#define _PGS_CONSTRAINTMANAGER_H_

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
      RefMat getViewLin(RefMat J, size_t i);
      //const Eigen::ConstRowsBlockXpr getViewLin(const ConstRefMat J, const Index i) const;
      RefMat getViewNonLin(RefMat J, size_t i);
      //const Eigen::ConstRowsBlockXpr getViewNonLin(const ConstRefMat J, const Index i) const;
    private:
      std::vector<Index> startLin;
      std::vector<Index> dimLin;
      std::vector<Index> startNonLin;
      std::vector<Index> dimNonLin;
  };
}


#endif // _PGS_CONSTRAINTMANAGER_H_
