#include <iostream>
#include <stdexcept>

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
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct SO3R2R3R2SO3(SO3R2R3R2, RotSpace);
  {
    std::cout << "Test RealSpace" << std::endl;
    Eigen::VectorXd v(3);
    v << 1,2,3;
    Point x = R3.createPoint(v);
    //std::cout << x << std::endl;
    //Eigen::MatrixXd Jac = R3.diffMap(x.value());
    //std::cout << Jac << std::endl;
    //Eigen::MatrixXd JacF(1,3);
    //JacF << 4,5,6;
    //R3.applyDiffMap(JacF, x.value());
    //std::cout << "JacF = " << JacF << std::endl;
    Eigen::Vector3d vOut;
    R3.invMap(vOut, x.value());
    std::cout << "v = " << v << std::endl;
    std::cout << "R3.invMap(vOut, x) = "<< vOut << std::endl;
  }
  {
    SO3<ExpMapMatrix> RotSpace;
    Eigen::Vector3d v ;//= Eigen::Vector3d::Random();
    v << 0.680375, -0.211234,  0.566198;
    Point x = RotSpace.getIdentity();
    Point id = RotSpace.getIdentity();
    x = x+v;
    Eigen::Vector3d vOut;
    std::cout << "v = " << v << std::endl;
    RotSpace.invMap(vOut, x.value());
    std::cout << "RotSpace.invMap(vOut, x) = "<< vOut << std::endl;
    RotSpace.minus(vOut, x.value(), id.value());
    std::cout << "RotSpace.minus(vOut, x.value(), id.value()) = "<< vOut << std::endl;
    //std::cout << "v = "<< v.transpose() << std::endl;
    //std::cout << "x = "<< x << std::endl;
    //Eigen::MatrixXd Jinv = RotSpace.diffInvMap(x.value());
    //std::cout << "Jinv = " << std::endl << Jinv << std::endl;
    //Eigen::MatrixXd J = RotSpace.diffMap(x.value());
    //std::cout << "J = " << std::endl << J << std::endl;
    //Eigen::Matrix<double, 2, 3> out;
    //Eigen::Matrix<double, 2, 9> in = Eigen::Matrix<double, 2, 9>::Random();
    //in << 1,2,3,4,5,6,7,8,9,
    //9,8,7,6,5,4,3,2,1;
    //RotSpace.applyDiffMap(out,in,x.value());
    //std::cout << "in = "  << std::endl << in << std::endl;
    //std::cout << "out = "  << std::endl << out << std::endl;

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

