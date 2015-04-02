#include <iostream>

#include <manifolds/RealSpace.h>
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE PGSolver
#endif

#include <boost/test/unit_test.hpp>

using namespace pgs;

BOOST_AUTO_TEST_CASE(typeIdEqualTest)
{
  RealSpace r5(5);
  RealSpace r6(6);
  S2 s1;
  S2 s1_1;
  SO3<ExpMapMatrix> so3M;
  SO3<ExpMapQuaternion> so3Q;
  SO3<ExpMapMatrix> so3M_1;
  SO3<ExpMapQuaternion> so3Q_1;

  BOOST_CHECK_EQUAL(r5.getTypeId(), r6.getTypeId());
  BOOST_CHECK_EQUAL(s1.getTypeId(), s1_1.getTypeId());
  BOOST_CHECK_EQUAL(so3M.getTypeId(), so3M_1.getTypeId());
  BOOST_CHECK_EQUAL(so3Q.getTypeId(), so3Q_1.getTypeId());

}

BOOST_AUTO_TEST_CASE(typeIdDifferentTest)
{
  std::vector<pgs::Manifold*> manifolds;

  manifolds.push_back(new RealSpace(5));
  manifolds.push_back(new S2());
  manifolds.push_back(new SO3<ExpMapMatrix>(5));
  manifolds.push_back(new SO3<ExpMapQuaternion>(5));

  for (size_t i = 0; i < manifolds.size(); ++i)
    {
      for (size_t j = i+1; j < manifolds.size(); ++j)
	{
	  BOOST_CHECK(manifolds[i]->getTypeId() != manifolds[j]->getTypeId());
	}
    }

  for (size_t i = 0; i < manifolds.size(); ++i)
    {
      delete manifolds[i];
    }
}
