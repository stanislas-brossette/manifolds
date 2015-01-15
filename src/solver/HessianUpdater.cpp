#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/Manifold.h>
#include <pgsolver/manifolds/Point.h>

#include <pgsolver/solver/SolverOptions.h>
#include <pgsolver/solver/HessianUpdater.h>


namespace pgs
{
  void HessianUpdater::hessianUpdate(Eigen::MatrixXd& H, const Point& x, 
      const double alpha, const Eigen::VectorXd& step, 
      const Eigen::MatrixXd& prevDiffLag, const Eigen::MatrixXd& diffLag, 
      const SolverOptions& solverOptions)
  {
    Index dim = x.getManifold().dim();
    Eigen::VectorXd sk(dim);
    Eigen::VectorXd yk(dim);
    x.getManifold().applyTransport(sk, alpha*step, x.value(), alpha*step); 
    x.getManifold().applyTransport(yk, prevDiffLag.transpose(), x.value(), alpha*step);
    yk = diffLag.transpose() - yk;
    x.getManifold().applyTransport(H, H, x.value(), alpha*step);
    x.getManifold().applyInvTransport(H, H, x.value(), alpha*step);

    if(sk.transpose()*yk<=0)
    {
      std::cerr << "Warning: Secant equation not respected... <sk,yk> <= 0 " << std::endl;
    }

    switch(solverOptions.hessianUpdateMethod){
      case BFGS: computeBFGS(H, sk, yk); break;
      case SR1: computeSR1(H, sk, yk); break;
      case EXACT: 
        throw std::runtime_error("EXACT update is not implemented yet"); 
        break;
    }
  }

  void HessianUpdater::computeBFGS(Eigen::MatrixXd& B, const Eigen::VectorXd& s,const Eigen::VectorXd& y)
  {
    Eigen::VectorXd Bs;
    Bs.noalias() = B*s;
    double sBs = s.transpose()*B*s;
    double theta;
    if(s.transpose()*y >= 0.2*sBs)
      theta = 1;
    else
      theta = (0.8*sBs)/(sBs-s.transpose()*y);

    Eigen::VectorXd r = theta*y + (1-theta)*Bs;

    B = B - (Bs*Bs.transpose())/sBs + (r*r.transpose())/(s.transpose()*r);

    //TODO Implement EigenValue control (LDL)
  }

  void HessianUpdater::computeSR1(Eigen::MatrixXd& B, const Eigen::VectorXd& s,const Eigen::VectorXd& y)
  {
    double r = 1e-8;
    Eigen::VectorXd Bs;
    Bs.noalias() = B*s;
    Eigen::VectorXd ymBs = y - Bs;
    if(fabs(s.transpose()*ymBs)>=r*s.lpNorm<1>()*ymBs.lpNorm<1>())
      B = B + (ymBs*ymBs.transpose())/(ymBs.transpose()*s);
    //else B = B
  }

}
