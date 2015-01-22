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
      const SolverOptions& opt)
  {
    if(step.lpNorm<Eigen::Infinity>() == 0)
    {
      //Should be used on first iteration
      //Nothing to do here
      return;
    }
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

    switch(opt.hessianUpdateMethod){
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

  void HessianUpdater::hessianUpdateIndividually(
          Eigen::MatrixXd& H, Eigen::MatrixXd& HCost, std::vector<Eigen::MatrixXd>& HCstr,
          const Eigen::VectorXd& lagMultNonLinCstr,
          const Point& x, const double alpha, const Eigen::VectorXd& step,
          const Eigen::MatrixXd& prevDiffObj, const Eigen::MatrixXd& diffObj,
          const Eigen::MatrixXd& prevDiffCstr, const Eigen::MatrixXd& diffCstr,
          const SolverOptions& opt)
  {
    if(step.lpNorm<Eigen::Infinity>() == 0)
    {
      //Should be used on first iteration
      //Nothing to do here
      return;
    }
     
    Index dim = x.getManifold().dim();
    Eigen::VectorXd sk(dim);
    Eigen::VectorXd yk(dim);
    x.getManifold().applyTransport(sk, alpha*step, x.value(), alpha*step);

    //Update of the Cost Hessian
    x.getManifold().applyTransport(yk, prevDiffObj.transpose(), x.value(), alpha*step);
    yk = diffObj.transpose() - yk;
    if(sk.transpose()*yk<=0)
    {
      std::cerr << "Warning: Secant equation not respected... <sk,yk> <= 0 " << std::endl;
    }
    //std::cout << "sk = \n" << sk.transpose() << std::endl;
    //std::cout << "yk = \n" << yk.transpose() << std::endl;
    x.getManifold().applyTransport(HCost, HCost, x.value(), alpha*step);
    x.getManifold().applyInvTransport(HCost, HCost, x.value(), alpha*step);
    switch(opt.hessianUpdateMethod){
      case BFGS: computeBFGS(HCost, sk, yk); break;
      case SR1: computeSR1(HCost, sk, yk); break;
      case EXACT:
        throw std::runtime_error("EXACT update is not implemented yet");
        break;
    }
    //std::cout << "HCost = \n" << HCost << std::endl;

    //Update of the Cstr Hessians
    for(size_t i = 0; i < HCstr.size(); ++i)
    {
      x.getManifold().applyTransport(yk, prevDiffCstr.row((int)i).transpose(), x.value(), alpha*step);
      yk = diffCstr.row((int)i).transpose() - yk;
      if(sk.transpose()*yk<=0)
      {
        std::cerr << "Warning: Secant equation not respected... <sk,yk> <= 0 " << std::endl;
      }
      //std::cout << "sk = \n" << sk.transpose() << std::endl;
      //std::cout << "yk = \n" << yk.transpose() << std::endl;
      x.getManifold().applyTransport(HCstr[i], HCstr[i], x.value(), alpha*step);
      x.getManifold().applyInvTransport(HCstr[i], HCstr[i], x.value(), alpha*step);
      switch(opt.hessianUpdateMethod){
        case BFGS: computeBFGS(HCstr[i], sk, yk); break;
        case SR1: computeSR1(HCstr[i], sk, yk); break;
        case EXACT:
          throw std::runtime_error("EXACT update is not implemented yet");
          break;
      }
      //std::cout << "HCstr[i] = \n" << HCstr[i] << std::endl;
    }

    //Update of the global Hessian
    H = HCost;
    for(size_t i = 0; i < HCstr.size(); ++i)
    {
      H += lagMultNonLinCstr[(Index)i]*HCstr[i];
    }
    //std::cout << "H = \n" << H << std::endl;
  }

}
