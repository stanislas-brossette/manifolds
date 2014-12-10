#include <iostream>
#include <stdexcept>

#include <pgsolver/SO3.h>
#include <pgsolver/RealSpace.h>
#include <pgsolver/CartesianProduct.h>
#include <pgsolver/Point.h>
#include <pgsolver/ExpMapMatrix.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE PGSolver 
#endif

#include <boost/test/unit_test.hpp>

using namespace pgs;

BOOST_AUTO_TEST_CASE(CartProdConstructor)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(P, RotSpace);
  BOOST_CHECK_EQUAL(Q.dim(), 10);
  BOOST_CHECK_EQUAL(Q.representationDim(), 16);
  BOOST_CHECK_EQUAL(Q.numberOfSubmanifolds(), 2);
  BOOST_CHECK_EQUAL(P.numberOfSubmanifolds(), 3);
}

BOOST_AUTO_TEST_CASE(CartProdIdentity)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(P, RotSpace);
  Point x = Q.getIdentity();
  BOOST_CHECK_EQUAL(x.value().size(), 16);
  BOOST_CHECK_EQUAL(x.value()[0], 0);
  BOOST_CHECK_EQUAL(x.value()[1], 0);
  BOOST_CHECK_EQUAL(x.value()[2], 0);
  BOOST_CHECK_EQUAL(x.value()[3], 0);
  BOOST_CHECK_EQUAL(x.value()[4], 0);
  BOOST_CHECK_EQUAL(x.value()[5], 0);
  BOOST_CHECK_EQUAL(x.value()[6], 0);
  BOOST_CHECK_EQUAL(x.value()[7], 1);
  BOOST_CHECK_EQUAL(x.value()[8], 0);
  BOOST_CHECK_EQUAL(x.value()[9], 0);
  BOOST_CHECK_EQUAL(x.value()[10], 0);
  BOOST_CHECK_EQUAL(x.value()[11], 1);
  BOOST_CHECK_EQUAL(x.value()[12], 0);
  BOOST_CHECK_EQUAL(x.value()[13], 0);
  BOOST_CHECK_EQUAL(x.value()[14], 0);
  BOOST_CHECK_EQUAL(x.value()[15], 1);
}

BOOST_AUTO_TEST_CASE(CartProdIncrement)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(P, RotSpace);
  Point x = Q.getIdentity();
  Eigen::VectorXd vy(10);
  vy << 1,2,3,4,5,6,7,0.1,0.2,0.3;
  x.increment(vy);
  x.increment(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 16);
  BOOST_CHECK_EQUAL(x.value()[0], 2);
  BOOST_CHECK_EQUAL(x.value()[1], 4);
  BOOST_CHECK_EQUAL(x.value()[2], 6);
  BOOST_CHECK_EQUAL(x.value()[3], 8);
  BOOST_CHECK_EQUAL(x.value()[4], 10);
  BOOST_CHECK_EQUAL(x.value()[5], 12);
  BOOST_CHECK_EQUAL(x.value()[6], 14);
  BOOST_CHECK_CLOSE(x.value()[7],   0.751909095300295, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[8],   0.583715086608147, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[9],  -0.306446422838863, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[10], -0.507379423623623, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[11],  0.809160842538688, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[12],  0.296352579515415, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[13],  0.420949917315650, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[14], -0.067345590561841, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[15],  0.904580421269344, 1e-12);
}

BOOST_AUTO_TEST_CASE(CartProdAddition)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(RotSpace, P);
  Point x = Q.getIdentity();
  Eigen::VectorXd vy(10);
  vy << 0.1,0.2,0.3,1,2,3,4,5,6,7;
  x = x + vy + vy;
  BOOST_CHECK_EQUAL(x.value().size(), 16);
  BOOST_CHECK_CLOSE(x.value()[0],  0.751909095300295, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[1],  0.583715086608147, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[2], -0.306446422838863, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[3], -0.507379423623623, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[4],  0.809160842538688, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[5],  0.296352579515415, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[6],  0.420949917315650, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[7], -0.067345590561841, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[8],  0.904580421269344, 1e-12);
  BOOST_CHECK_EQUAL(x.value()[9], 2); 
  BOOST_CHECK_EQUAL(x.value()[10],4); 
  BOOST_CHECK_EQUAL(x.value()[11],6); 
  BOOST_CHECK_EQUAL(x.value()[12],8); 
  BOOST_CHECK_EQUAL(x.value()[13],10);
  BOOST_CHECK_EQUAL(x.value()[14],12);
  BOOST_CHECK_EQUAL(x.value()[15],14);
}

BOOST_AUTO_TEST_CASE(CartProSubstraction)
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  Point x = SO3R2R3R2.getIdentity();
  Point y = SO3R2R3R2.getIdentity();
  Eigen::VectorXd vx(10);
  Eigen::VectorXd vy(10);
  vx << 1,0.1,1,1,2,3,4,5,6,7;
  vy << 0.07,3,0.01,4,6,2,1,4,6,4;
  x = x + vx;
  y = y + vy;
  Point z = x+(y-x);

  BOOST_CHECK_EQUAL(z.value().size(), 16);
  for (int i = 0; i < 10; ++i)
  {
    BOOST_CHECK_CLOSE(y.value()[i], z.value()[i], 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(CardProdPointInvMap)
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct Space(SO3R2R3R2, RotSpace);
  Point x = Space.getIdentity();
  Eigen::VectorXd vy = Eigen::VectorXd::Random(Space.dim());;
  x = x + vy;
  Eigen::VectorXd z(Space.dim());
  Space.invMap(z, x.value()); 
  BOOST_CHECK(z.isApprox(vy));
}

BOOST_AUTO_TEST_CASE(CartProInvMap)
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2SO3(R2, RotSpace);
  CartesianProduct R2R3(R2, R3);
  CartesianProduct R2SO3R2R3(R2SO3, R2R3);
  Point x = R2SO3R2R3.getIdentity();
  Point Id = R2SO3R2R3.getIdentity();
  Eigen::VectorXd vx(10);
  Eigen::VectorXd vy(10);
  vx << -7,2,1,0.1,1,3,4,5,6,7;
  vy << 4,6,0.07,3,0.01,2,1,4,6,4;
  x = x + vx +vy;
  Eigen::VectorXd x0 = x.invMap();

  Point newX = Id.increment(x0); 

  BOOST_CHECK_EQUAL(newX.value().size(), 16);
  for (int i = 0; i < 10; ++i)
  {
    BOOST_CHECK_CLOSE(newX.value()[i], x.value()[i], 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(CardProdDiff)
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2SO3(R2, RotSpace);
  Eigen::MatrixXd J;
  Eigen::MatrixXd Jtest(11,5);
  Jtest <<  1, 0, 0                 , 0                 , 0                   , 
            0, 1, 0                 , 0                 , 0                   , 
            0, 0, 0                 , -0.573417053935868, -0.109861317404411  ,
            0, 0, 0                 ,  0.249497069246539,  0.920480138494529  ,
            0, 0, 0                 , -0.780348700705586,  0.375029072973363  ,
            0, 0,  0.573417053935868, 0                 ,  -0.811864134688605 ,
            0, 0, -0.249497069246539, 0                 ,  -0.300778202459022 ,
            0, 0,  0.780348700705586, 0                 ,   0.500408932506048 ,
            0, 0,  0.109861317404411,  0.811864134688605, 0                   ,
            0, 0, -0.920480138494529,  0.300778202459022, 0                   ,
            0, 0, -0.375029072973363, -0.500408932506048, 0                   ;
  Point x = R2SO3.getIdentity();
  Eigen::VectorXd v(5);
  v << 243.27432598, -2314327.23748, 0.3403857, 0.58526775, 0.223811;
  x.increment(v);
  J = R2SO3.diffMap(x.value());
  BOOST_CHECK(J.isApprox(Jtest));
}

BOOST_AUTO_TEST_CASE(CardProdApplyDiff)
{
  int c = 5;
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct Space(SO3R2R3R2, RotSpace);
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
  Point x = Space.getIdentity();
  x.increment(Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*Space.diffMap(x.value());
  Eigen::MatrixXd J(c,dim);
  Space.applyDiffMap(J, Jf, x.value());
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(CardProdDiffInv)
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2SO3(R2, RotSpace);
  CartesianProduct R2SO3R3(R2SO3, R3);
  Eigen::MatrixXd J;
  Eigen::MatrixXd Jtest(8,14);
  Jtest << 1,0,0                 ,0                ,0                 ,                 0,                  0,                 0,                 0,                  0,                  0,0,0,0,
           0,1,0                 ,0                ,0                 ,                 0,                  0,                 0,                 0,                  0,                  0,0,0,0,
           0,0,-0.064043491813865,                0,                 0,                 0, -0.064043491813865, 0.545030346992499,                 0, -0.545030346992499, -0.064043491813865,0,0,0,
           0,0,-0.110117993664377,                0,-0.545030346992499,                 0, -0.110117993664377,                 0, 0.545030346992499,                  0, -0.110117993664377,0,0,0,
           0,0,-0.042109988599266,0.545030346992499,                 0,-0.545030346992499, -0.042109988599266,                 0,                 0,                  0, -0.042109988599266,0,0,0,
           0,0,                 0,                0,                 0,                 0,                  0,                 0,                 0,                  0,                  0,1,0,0, 
           0,0,                 0,                0,                 0,                 0,                  0,                 0,                 0,                  0,                  0,0,1,0, 
           0,0,                 0,                0,                 0,                 0,                  0,                 0,                 0,                  0,                  0,0,0,1; 
  Point x = R2SO3R3.getIdentity();
  Eigen::VectorXd v(8);
  v << 243.27432598, -2314327.23748, 0.3403857, 0.58526775, 0.223811, 3.08, 0.00000001, 232.5;
  x.increment(v);
  J = R2SO3R3.diffInvMap(x.value());
  BOOST_CHECK(J.isApprox(Jtest));
}

BOOST_AUTO_TEST_CASE(CardProdApplyInvDiff)
{
  int c = 5;
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct Space(SO3R2R3R2, RotSpace);
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,dim);
  Point x = Space.getIdentity();
  x.increment(Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*Space.diffInvMap(x.value());
  Eigen::MatrixXd J(c,repDim);
  Space.applyDiffInvMap(J, Jf, x.value());
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(CardProdTransport)
{
  std::cout << "CardProdTransport" << std::endl;
  int c = 4;
  SO3<ExpMapMatrix> S;
  RealSpace R3(3);
  RealSpace R2(2);
  CartesianProduct Space(R3, S);
  Space.multiply(R2);
  Index dim = Space.dim();
  std::cout << "dim = " << dim<< std::endl;
  Eigen::MatrixXd H(dim,c);
  H <<   1, 2, 3, 4,
         5, 6, 7, 8,
         9,10,11,12,
        13,14,15,16,
        17,18,19,20,
        21,22,23,24,
        25,26,27,28,
        29,30,31,32;
  std::cout << "H = " << H << std::endl;
  Eigen::MatrixXd Hout(dim,c);
  Eigen::VectorXd v(dim);
  v <<  0.141886338627215,
        0.421761282626275,
        0.915735525189067,
        0.287150084472884,
        0.145612694616852,
        0.240084140666640,
        0.792207329559554,
        0.959492426392903;
  Point x = Space.getIdentity();
  x.increment(Eigen::VectorXd::Random(8));
  Eigen::MatrixXd expectedRes(dim,c);
  expectedRes <<   1, 2, 3, 4,
                   5, 6, 7, 8,
                   9,10,11,12,
                  12.562951305087033, 13.486740548139206, 14.410529791191379, 15.334319034243553,
                  13.623940511618414, 14.546891198380022, 15.469841885141630, 16.392792571903236,
                  23.570330914984936, 24.708212923027673, 25.846094931070407, 26.983976939113141,
                  25,26,27,28,
                  29,30,31,32;
                 
  std::cout << "expectedRes = " << expectedRes << std::endl;
  Space.applyTransport(Hout, H, x.value(), v);
  std::cout << "Hout = " << Hout << std::endl;
  BOOST_CHECK(expectedRes.isApprox(Hout));
}

BOOST_AUTO_TEST_CASE(CardProdInvTransport)
{
  std::cout << "SO3InvTransport" << std::endl;
  int r = 4;
  SO3<ExpMapMatrix> S;
  RealSpace R3(3);
  RealSpace R2(2);
  CartesianProduct Space(R3, S);
  Space.multiply(R2);
  Index dim = Space.dim();
  Eigen::MatrixXd H(r,dim);
  H <<  1, 2, 3, 4, 5, 6, 7, 8,
        9,10,11,12,13,14,15,16,
       17,18,19,20,21,22,23,24,
       25,26,27,28,29,30,31,32;
  std::cout << "H = " << H << std::endl;
  Eigen::MatrixXd Hout(r,dim);
  Eigen::VectorXd v(dim);
  v <<  0.013851417189346,
        0.029139534370754,
        0.247037348498188,
        0.208448586892745,
        0.095129844018258,
        0.285066614651506,
        0.010333824150873,
        0.131623307896919;

  Point x = Space.getIdentity();
  x.increment(Eigen::VectorXd::Random(8));
  Eigen::MatrixXd expectedRes(r,dim);

  expectedRes <<   1, 2, 3, 3.211060456124126,  4.703367298475802,  6.675883971635843, 7, 8,
                   9,10,11, 9.681461929634750, 12.995126625170094, 15.697005411887588,15,16,
                  17,18,19,16.151863403145374, 21.286885951864392, 24.718126852139335,23,24,
                  25,26,27,22.622264876655997, 29.578645278558682, 33.739248292391082,31,32;
  std::cout << "expectedRes = " << expectedRes << std::endl;
  Space.applyInvTransport(Hout, H, x.value(), v);
  std::cout << "Hout = " << Hout << std::endl;
  BOOST_CHECK(expectedRes.isApprox(Hout));
}

#if   EIGEN_WORLD_VERSION > 3 \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && EIGEN_MINOR_VERSION > 0)
BOOST_AUTO_TEST_CASE(CardProdNoAllocation)
{
  //We only test here that the operations on the manifold do not create
  //temporary. Passing arguments that are not recognize by the Eigen::Ref will
  //create temporaries, but this is the user's fault.
  const int r = 100;
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct S(SO3R2R3R2, RotSpace);
  Index dim = S.dim();
  Index repDim = S.representationDim();
  Eigen::VectorXd x = Eigen::VectorXd::Random(repDim);
  Eigen::VectorXd p = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd y = Eigen::VectorXd::Random(repDim);
  Eigen::VectorXd z(repDim);
  Eigen::VectorXd d(dim);
  Eigen::MatrixXd J0 = Eigen::MatrixXd::Random(r, repDim);
  Eigen::MatrixXd J1(r, dim);
  Eigen::MatrixXd J2(r, repDim);
  Eigen::MatrixXd H0 = Eigen::MatrixXd::Random(S.dim(), S.dim());
  Eigen::MatrixXd H1 = Eigen::MatrixXd::Random(S.dim(), S.dim());
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Random(S.dim(), S.dim());

  //The first call to the following methods might trigger a memory allocation
  //depending on the size of the Ji and the initial buffer size inside S.
  //However, subsequent calls should not require any allocation, what we check
  //after.
  S.applyDiffMap(J1, J0, x);
  S.applyDiffInvMap(J2, J1, x);

  Eigen::internal::set_is_malloc_allowed(false);
  {
    S.plus(z, x, p);
    S.minus(d, x, y);
    S.invMap(d, x);
    S.applyDiffMap(J1, J0, x);
    S.applyDiffInvMap(J2, J1, x);
    S.applyTransport(H1, H0, x, p);
    S.applyInvTransport(H2, H0, x, p);
  }
  Eigen::internal::set_is_malloc_allowed(true);
}
#endif
