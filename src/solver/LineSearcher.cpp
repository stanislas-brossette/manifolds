#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/Manifold.h>
#include <pgsolver/manifolds/Point.h>

#include <pgsolver/solver/SolverOptions.h>
#include <pgsolver/solver/LineSearcher.h>

namespace pgs
{
  double LineSearcher::LineSearch(Solver& s, 
      Problem& p, ProblemEvaluation& probE,  
      Filter& filter_, Eigen::VectorXd& step, SolverOptions& opt,
      Eigen::VectorXi feasibilityStatus)
  {
    if(opt.VERBOSE >= 2) 
    {
      std::cout << "---- Line-search ----"<< std::endl;
      std::cout << "Starting point: "<< p.x() << std::endl;
      std::cout << "Search direction: "<< step.transpose() << std::endl;
      std::cout << "Filter = " << std::endl;
      filter_.print();
    }

    double alpha = 1.0;

    Eigen::Vector2d FH;

    while(alpha > 1e-4)
    {
      p.setZ(alpha*step);
      s.updateAllCstr(p);

      FH = s.computeFHforFilter(probE, feasibilityStatus);

      if(opt.VERBOSE >= 2) 
      {
        std::cout << "Trial with alpha = " << alpha << std::endl;
        std::cout << "F = " << FH[0] << std::endl;
        std::cout << "H = " << FH[1] << std::endl;
      }

      if(filter_.accepts(FH))
      {
        filter_.add(FH);
        if(opt.VERBOSE >= 2) 
        {
          std::cout << "Filter accepted: " << FH.transpose() << std::endl;
          std::cout << "alpha = " << alpha << std::endl;
          std::cout << "---------------------"<< std::endl;
        }
        return alpha;
      }
      else
      {
        if(opt.VERBOSE >= 2) 
        {
          std::cout << "Filter refused: " << FH.transpose() << std::endl;
        }
        alpha = alpha*0.9;
      }
    }
    std::cout << "\033[1;31mWARNING: Line-Search did not converge\033[0m" << std::endl;
    //throw std::runtime_error("Line-Search did not converge");
    std::cout << "---------------------"<< std::endl;
    return alpha;
  }
}
