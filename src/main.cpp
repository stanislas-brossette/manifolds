#include <iostream>
#include <stdexcept>
#define _USE_MATH_DEFINES
#include <math.h>

#include <Eigen/Core>

#include <pgsolver/pgs_assert.h>
#include <pgsolver/RealSpace.h>
#include <pgsolver/SO3.h>
#include <pgsolver/CartesianProduct.h>
#include <pgsolver/Point.h>
#include <pgsolver/ExpMapMatrix.h>
#include <pgsolver/ReusableTemporaryMap.h>

using namespace pgs;

int main()
{
  std::cout << "Using: Eigen" << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION <<"." << EIGEN_MINOR_VERSION<< std::endl;
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct SO3R2R3R2SO3(SO3R2R3R2, RotSpace);
  CartesianProduct R3SO3R3SO3(R3, RotSpace);
  R3SO3R3SO3.multiply(R3);
  R3SO3R3SO3.multiply(RotSpace);
  {
    std::cout << "Test RealSpace" << std::endl;
    Eigen::VectorXd v(3);
    Eigen::MatrixXd M(3,4);
    M = Eigen::MatrixXd::Random(3,4);
    Eigen::MatrixXd MTrans(3,4);
    v << 1,2,3;
    Point x = R3.createPoint(v);
    R3.applyTransport(MTrans, M, v);
    std::cout << "v = " << std::endl<< v << std::endl;
    std::cout << "M = " << std::endl<< M << std::endl;
    std::cout << "R3.applyTransport(MTrans, M, v) = "<< std::endl<< MTrans << std::endl;
  }
  {
    std::cout << "Test RotSpace" << std::endl;
    Eigen::VectorXd v(3);
    Eigen::MatrixXd M(3,4);
    M = Eigen::MatrixXd::Random(3,4);
    Eigen::MatrixXd MTrans(3,4);
    v << 0,0,0;
    Point x = RotSpace.getIdentity();
    x.increment(v);
    RotSpace.applyTransport(MTrans, M, v);
    std::cout << "v = " << std::endl<< v << std::endl;
    std::cout << "M = " << std::endl<< M << std::endl;
    std::cout << "RotSpace.applyTransport(MTrans, M, v) = " << std::endl << MTrans << std::endl;
  }
  {
    std::cout << "Test CardProdSpace" << std::endl;
    Eigen::VectorXd v(12);
    Eigen::MatrixXd M(12,12);
    //M = Eigen::MatrixXd::Random(12,4);
    for(Index i = 0; i<12; ++i)
    {
      for(Index j = 0; j<12; ++j)
      {
        std::cout << 100*j+i << std::endl;
        M(i,j) = 100*(double)j+(double)i;
      }
    }
    Eigen::MatrixXd MTrans(12,12);
    Eigen::MatrixXd MInvTrans(12,12);
    Eigen::MatrixXd MDoubleTrans(12,12);
    v << 1,2,3,0,M_PI/2,0,4,5,6,M_PI/2,0,0;
    Point x = R3SO3R3SO3.getIdentity();
    x.increment(v);
    R3SO3R3SO3.applyTransport(MTrans, M, v);
    R3SO3R3SO3.applyInvTransport(MInvTrans, M, v);
    R3SO3R3SO3.applyInvTransport(MDoubleTrans, MTrans, v);
    std::cout << "v = " << std::endl << v << std::endl;
    std::cout << "M = "<< std::endl << M << std::endl;
    std::cout << "MTrans = "<< std::endl << MTrans << std::endl;
    std::cout << "MInvTrans = "<< std::endl << MInvTrans << std::endl;
    std::cout << "MDoubleTrans = "<< std::endl << MDoubleTrans << std::endl;
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

  //ReusableTemporaryMap rtm;
  //Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> tmp = rtm.getMap(5, 7);
  //tmp.setIdentity();
  //std::cout << "initial 5x7 matrix: " << std::endl << tmp << std::endl;
  //Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> tmp2 = rtm.getMap(4, 6);
  //std::cout << "new map 4x6, smaller. No reallocation, only use initialized values of previous map : " << std::endl << tmp2 << std::endl;
  //Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> tmp3 = rtm.getMap(7, 7);
  //std::cout << "bigger 7x7, no reallocation but some values are not initialized : " << std::endl << tmp3 << std::endl;
  //Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> tmp4 = rtm.getMap(17, 17);
  //std::cout << "a 17x17 matrix doesn't fit in the initial buffer. Memory is reallocated and initialized values may be lost : " << std::endl << tmp4 << std::endl;
#ifdef _WIN32
  system("pause");
#endif
  return 0;
}

