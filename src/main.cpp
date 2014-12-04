#include <iostream>
#include <stdexcept>

#include <Eigen/Core>

#include <pgsolver/pgs_assert.h>
#include <pgsolver/RealSpace.h>
#include <pgsolver/SO3.h>
#include <pgsolver/CartesianProduct.h>
#include <pgsolver/Point.h>
#include <pgsolver/ExpMapMatrix.h>

using namespace pgs;

int main()
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct SO3R2R3R2SO3(SO3R2R3R2, RotSpace);
  //{
  //  std::cout << "Test RealSpace" << std::endl;
  //  Eigen::VectorXd v(3);
  //  v << 1,2,3;
  //  Point x = R3.createPoint(v);
  //  std::cout << x << std::endl;
  //  Eigen::MatrixXd Jac = R3.diffMap(x.value());
  //  std::cout << Jac << std::endl;
  //  Eigen::MatrixXd JacF(1,3);
  //  JacF << 4,5,6;
  //  R3.applyDiffMap(JacF, x.value());
  //  std::cout << "JacF = " << JacF << std::endl;
  //}
  {
    Index c = 5;
    SO3<ExpMapMatrix> RotSpace;
    Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,9);
    Eigen::MatrixXd Jres = Eigen::MatrixXd::Random(c,3);
    Point x = RotSpace.getIdentity();
    RotSpace.applyDiffMap(Jres, Jf, x.value());
    std::cout << "Jf=" << std::endl << Jf << std::endl;
    std::cout << "Jres=" << std::endl << Jres << std::endl;
    
    bool worked = true;

    for (int i = 0; i<13; ++i)
    {
      Eigen::MatrixXd G = Eigen::MatrixXd::Zero(c,15);
      G.middleCols(3,9) = Jf;
      Eigen::Map<Eigen::MatrixXd> Gf(G.data()+3*c,c,9);
      Eigen::Map<Eigen::MatrixXd> Gres(G.data()+i*c,c,3);
      try
      {
        RotSpace.applyDiffMap(Gres,Gf,x.value());
      }
      catch (pgs_exception &e)
      {
        std::cout << "At iteration " << i << ", exception caught" << std::endl;
        continue;
      }
      worked = worked && Jres.isApprox(Gres);
      std::cout << "At iteration " << i << ": Gres-Jres=" << std::endl << Gres-Jres << std::endl;
    }
  }
  
  //{
  //  std::cout << "Test SO3" << std::endl;
  //  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(5,9);
  //  Point x = RotSpace.getIdentity();
  //  Eigen::MatrixXd Jac = RotSpace.diffMap(x.value());
  //  Eigen::MatrixXd expectedRes;
  //  expectedRes = Jf*Jac;
  //  std::cout << "Jf=" << std::endl << Jf << std::endl;
  //  Eigen::Map<Eigen::MatrixXd> J(Jf.data(),5,3);
  //  RotSpace.applyDiffMap(J, Jf, x.value());
  //  std::cout << "J after applyDiffMap = " << std::endl << J << std::endl;
  //  std::cout << "Jf*Jac - J=" << std::endl << expectedRes - J << std::endl;
  //}
  //{
  //  std::cout << "Test CardProd" << std::endl;
  //  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(5,25);
  //  Point x = SO3R2R3R2SO3.getIdentity();
  //  Eigen::MatrixXd expectedRes = Jf*SO3R2R3R2SO3.diffMap(x.value());
  //  std::cout << "Jf=" << std::endl << Jf << std::endl;
  //  Eigen::Map<Eigen::MatrixXd> J(Jf.data(),5,13);
  //  SO3R2R3R2SO3.applyDiffMap(J, Jf, x.value());
  //  std::cout << "J after applyDiffMap = " << std::endl << J << std::endl;
  //  std::cout << "Jf*Jac - J=" << std::endl << expectedRes - J << std::endl;
  //}

#ifdef _WIN32
  system("pause");
#endif
  return 0;
}

