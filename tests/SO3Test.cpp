#include <iostream>
#include <stdexcept>
#include <pgsolver/pgs_assert.h>
#include <pgsolver/SO3.h>
#include <pgsolver/Point.h>
#include <pgsolver/ExpMapMatrix.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE PGSolver 
#endif

#include <boost/test/unit_test.hpp>

using namespace pgs;

BOOST_AUTO_TEST_CASE(RotSpaceConstructor)
{
  SO3<ExpMapMatrix> RotSpace;
  BOOST_CHECK_EQUAL(RotSpace.dim(), 3);
  BOOST_CHECK_EQUAL(RotSpace.representationDim(), 9);
  BOOST_CHECK_EQUAL(RotSpace.numberOfSubmanifolds(), 1);
}

BOOST_AUTO_TEST_CASE(SO3Constructor)
{
  SO3<ExpMapMatrix> RotSpace;
  Point x = RotSpace.getIdentity();
  Eigen::VectorXd v(9);
  v << 1,0,0,0,1,0,0,0,1;
  Point y = RotSpace.createPoint(v);
  BOOST_CHECK_EQUAL(x.value().size(), 9);
  BOOST_CHECK_EQUAL(x.value()[0], 1);
  BOOST_CHECK_EQUAL(x.value()[1], 0);
  BOOST_CHECK_EQUAL(x.value()[2], 0);
  BOOST_CHECK_EQUAL(x.value()[3], 0);
  BOOST_CHECK_EQUAL(x.value()[4], 1);
  BOOST_CHECK_EQUAL(x.value()[5], 0);
  BOOST_CHECK_EQUAL(x.value()[6], 0);
  BOOST_CHECK_EQUAL(x.value()[7], 0);
  BOOST_CHECK_EQUAL(x.value()[8], 1);
  BOOST_CHECK_EQUAL(y.value().size(), 9);
  BOOST_CHECK_EQUAL(y.value()[0], 1);
  BOOST_CHECK_EQUAL(y.value()[1], 0);
  BOOST_CHECK_EQUAL(y.value()[2], 0);
  BOOST_CHECK_EQUAL(y.value()[3], 0);
  BOOST_CHECK_EQUAL(y.value()[4], 1);
  BOOST_CHECK_EQUAL(y.value()[5], 0);
  BOOST_CHECK_EQUAL(y.value()[6], 0);
  BOOST_CHECK_EQUAL(y.value()[7], 0);
  BOOST_CHECK_EQUAL(y.value()[8], 1);
}

BOOST_AUTO_TEST_CASE(SO3Increment)
{
  SO3<ExpMapMatrix> RotSpace;
  Point x = RotSpace.getIdentity();
  Eigen::Vector3d vy;
  vy << 0.1,0.2,0.3;
  x.increment(vy);
  x.increment(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 9);
  BOOST_CHECK_CLOSE(x.value()[0],  0.751909095300295, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[1],  0.583715086608147, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[2], -0.306446422838863, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[3], -0.507379423623623, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[4],  0.809160842538688, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[5],  0.296352579515415, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[6],  0.420949917315650, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[7], -0.067345590561841, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[8],  0.904580421269344, 1e-12);
}

BOOST_AUTO_TEST_CASE(SO3Addition)
{
  SO3<ExpMapMatrix> RotSpace;
  Point y = RotSpace.getIdentity();
  Eigen::Vector3d vy;
  vy << 0.1,0.2,0.3;
  y=(y+vy)+vy;
  BOOST_CHECK_EQUAL(y.value().size(), 9);
  BOOST_CHECK_CLOSE(y.value()[0],  0.751909095300295, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[1],  0.583715086608147, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[2], -0.306446422838863, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[3], -0.507379423623623, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[4],  0.809160842538688, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[5],  0.296352579515415, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[6],  0.420949917315650, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[7], -0.067345590561841, 1e-12);
  BOOST_CHECK_CLOSE(y.value()[8],  0.904580421269344, 1e-12);
}

BOOST_AUTO_TEST_CASE(SO3Substraction)
{
  SO3<ExpMapMatrix> RotSpace;
  Eigen::Vector3d v( 0.2, 0.4, 0.6);
  Point R1 = RotSpace.getIdentity();
  R1 = R1 + v;
  Point R2 = R1 + v;
  Eigen::Vector3d d = R2-R1;
  BOOST_CHECK_CLOSE(d[0], 0.2, 1e-8);
  BOOST_CHECK_CLOSE(d[1], 0.4, 1e-8);
  BOOST_CHECK_CLOSE(d[2], 0.6, 1e-8);
}

BOOST_AUTO_TEST_CASE(SO3Diff)
{
  SO3<ExpMapMatrix> RotSpace;
  Eigen::MatrixXd J;
  Eigen::MatrixXd Jtest(9,3);
  Jtest << 0, 0, 0,  
           0, 0, 1,
           0,-1, 0,
           0, 0,-1,
           0, 0, 0,
           1, 0, 0,
           0, 1, 0,
          -1, 0, 0,
           0, 0, 0;
  Point x = RotSpace.createPoint();
  J = RotSpace.diffMap(x.value());
  BOOST_CHECK(J.isApprox(Jtest));
}

BOOST_AUTO_TEST_CASE(SO3ApplyDiff)
{
  SO3<ExpMapMatrix> RotSpace;
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(5,9);
  Point x = RotSpace.getIdentity();
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*RotSpace.diffMap(x.value());
  Eigen::MatrixXd J(5,3);
  RotSpace.applyDiffMap(J, Jf, x.value());
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(SO3ApplyDiffMemoryTest)
{
  Index c = 6;
  SO3<ExpMapMatrix> Space;
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
  Eigen::MatrixXd Jres = Eigen::MatrixXd::Random(c,dim);
  Point x = Space.getIdentity();
  Space.applyDiffMap(Jres, Jf, x.value());
  
  for (int i = 0; i<dim+repDim; ++i)
  {
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(c,repDim+2*dim);
    G.middleCols(dim,repDim) = Jf;
    Eigen::Map<Eigen::MatrixXd> Gf(G.data()+dim*c,c,repDim);
    Eigen::Map<Eigen::MatrixXd> Gres(G.data()+i*c,c,dim);
    try
    {
      Space.applyDiffMap(Gres,Gf,x.value());
    }
    catch (pgs_exception&)
    {
      ExpMapMatrix::applyDiffMapNoAssert_(Gres,Gf,x.value());
      BOOST_CHECK(!Jres.isApprox(Gres));
      continue;
    }
    BOOST_CHECK(Jres.isApprox(Gres));
  }
}

BOOST_AUTO_TEST_CASE(SO3ApplyDiffGuaranteedResultTest)
{
  Index c = 5;
  SO3<ExpMapMatrix> Space;
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
  Eigen::MatrixXd Jres = Eigen::MatrixXd::Random(c,dim);
  Point x = Space.getIdentity();
  Space.applyDiffMap(Jres, Jf, x.value());
  
  for (int i = 0; i<dim+1; ++i)
  {
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(c,repDim+dim);
    G.middleCols(dim,repDim) = Jf;
    Eigen::Map<Eigen::MatrixXd> Gf(G.data()+dim*c,c,repDim);
    Eigen::Map<Eigen::MatrixXd> Gres(G.data()+i*c,c,dim);
    Space.applyDiffMap(Gres,Gf,x.value());
    BOOST_CHECK(Jres.isApprox(Gres));
  }
}
