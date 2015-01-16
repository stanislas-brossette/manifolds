#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/Manifold.h>
#include <pgsolver/manifolds/Point.h>

#include <pgsolver/solver/SolverOptions.h>
#include <pgsolver/solver/LineSearcher.h>

namespace pgs
{
  double LineSearcher::LineSearch(Solver& s, Problem& p, Filter& filter_, Eigen::VectorXd& step)
  {
    double alpha = 1.0;
    std::cout << "In line-search"<< std::endl;
    std::cout << "Starting point: "<< p.x() << std::endl;
    std::cout << "Search direction: "<< step.transpose() << std::endl;
    std::cout << filter_.size() << std::endl;
    std::cout << "F = " << s.probEval().obj << std::endl;
    std::cout << "H = " << s.probEval().violCstr.transpose() << std::endl;

    Eigen::Vector2d FH;

    while(alpha > 1e-12)
    {
      p.setZ(alpha*step);
      s.updateObj(p);
      s.updateViolations(p);
      FH << s.probEval().obj, s.probEval().violCstr.lpNorm<1>();
      
      if(filter_.accepts(FH))
      {
        filter_.add(FH);
        return alpha;
      }
      else
      {
        alpha = alpha*0.9;
      }
    }
    std::cout << "\033[1;31mWARNING: Line-Search did not converge\033[0m" << std::endl;
    //throw std::runtime_error("Line-Search did not converge"); 
    return alpha;
  }
}
