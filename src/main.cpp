#include <iostream>
#include <stdexcept>
#define _USE_MATH_DEFINES
#include <math.h>

#include <Eigen/Core>

#include <pgsolver/manifolds/pgs_assert.h>
#include <pgsolver/manifolds/RealSpace.h>
#include <pgsolver/manifolds/SO3.h>
#include <pgsolver/manifolds/CartesianProduct.h>
#include <pgsolver/manifolds/Point.h>
#include <pgsolver/manifolds/ExpMapMatrix.h>
#include <pgsolver/manifolds/ReusableTemporaryMap.h>

#include <pgsolver/solver/ExampleGeometricProblem.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Solver.h>

#include <EigenQP/QLD.h>
#include <EigenQP/QuadProg.h>
#include <EigenQP/LSSOL.h>

using namespace pgs;

struct QP1
{
	QP1()
	{
		nrvar = 6;
		nreq = 3;
		nrineq = 2;

		Q.resize(nrvar, nrvar);
		Aeq.resize(nreq, nrvar);
		Aineq.resize(nrineq, nrvar);

		C.resize(nrvar);
		Beq.resize(nreq);
		Bineq.resize(nrineq);
		XL.resize(nrvar);
		XU.resize(nrvar);
		X.resize(nrvar);


		Aeq << 1., -1., 1., 0., 3., 1.,
					 -1., 0., -3., -4., 5., 6.,
					 2., 5., 3., 0., 1., 0.;
		Beq << 1., 2., 3.;

		Aineq << 0., 1., 0., 1., 2., -1.,
						 -1., 0., 2., 1., 1., 0.;
		Bineq << -1., 2.5;

		//with  x between ci and cs:
		XL << -1000., -10000., 0., -1000., -1000.,-1000.;
		XU << 10000., 100., 1.5, 100., 100., 1000.;

		//and minimize 0.5*x'*Q*x + p'*x with
		C << 1., 2., 3., 4., 5., 6.;
		Q.setIdentity();

		X << 1.7975426, -0.3381487, 0.1633880, -4.9884023, 0.6054943, -3.1155623;
	}

	int nrvar, nreq, nrineq;
	Eigen::MatrixXd Q, Aeq, Aineq;
	Eigen::VectorXd C, Beq, Bineq, XL, XU, X;
};

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
  CartesianProduct R3SO3(R3, RotSpace);
  CartesianProduct R3SO3R3SO3(R3, RotSpace);
  R3SO3R3SO3.multiply(R3);
  R3SO3R3SO3.multiply(RotSpace);
  {
    ExampleGeometricProblem myProb;
    Solver mySolver;
    Eigen::VectorXd v0(3);
    v0 << 1,2,3;
    Point x0 = myProb.M().createPoint(v0);
    mySolver.solve(myProb, x0);
  }
  {
    std::cout << "==================================================" << std::endl;
    std::cout << "================== Test LSSOL ====================" << std::endl;
    std::cout << "==================================================" << std::endl;
    double inf = std::numeric_limits<double>::infinity();

    QP1 qp1;

    int nrconstr = qp1.nreq + qp1.nrineq;
    Eigen::LSSOL lssol(qp1.nrvar, nrconstr);

    Eigen::VectorXd AL(nrconstr);
    AL.segment(0, qp1.nreq) = qp1.Beq;
    AL.segment(qp1.nreq, qp1.nrineq).fill(-inf);

    Eigen::VectorXd AU(nrconstr);
    AU.segment(0, qp1.nreq) = qp1.Beq;
    AU.segment(qp1.nreq, qp1.nrineq) = qp1.Bineq;

    Eigen::MatrixXd A(nrconstr, qp1.nrvar);
    A.block(0, 0, qp1.nreq, qp1.nrvar) = qp1.Aeq;
    A.block(qp1.nreq, 0, qp1.nrineq, qp1.nrvar) = qp1.Aineq;

    lssol.solve(qp1.Q, qp1.C,
      A, static_cast<int>(A.rows()), AL, AU,
      qp1.XL, qp1.XU);

    std::cout << lssol.result() - qp1.X << std::endl;
      
  }
#ifdef _WIN32
  system("pause");
#endif
  return 0;
}

