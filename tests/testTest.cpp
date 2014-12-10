#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include <pgsolver/defs.h>

using namespace pgs;

void foo1(RefVec out, RefVec in, RefVec p)
{
  out = in + p;
}
void foo2(RefVec out, const RefVec in, const RefVec p)
{
  out = in + p;
}
void foo3(RefVec out, ConstRefVec in, ConstRefVec p)
{
  out = in + p;
}
void foo4(RefVec out, const ConstRefVec in, const ConstRefVec p)
{
  out = in + p;
}

BOOST_AUTO_TEST_CASE(testAllocationFOO1)
{
  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::VectorXd y = Eigen::VectorXd::Random(4);
  Eigen::VectorXd z(4);
  Eigen::internal::set_is_malloc_allowed(false);
  {
    foo1(z, x, y);
  }
  Eigen::internal::set_is_malloc_allowed(true);
  std::cout << "x=" << x.transpose() << std::endl;
  std::cout << "y=" << y.transpose() << std::endl;
  std::cout << "z=" << z.transpose() << std::endl;
}
BOOST_AUTO_TEST_CASE(testAllocationFOO2)
{
  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::VectorXd y = Eigen::VectorXd::Random(4);
  Eigen::VectorXd z(4);
  Eigen::internal::set_is_malloc_allowed(false);
  {
    foo2(z, x, y);
  }
  Eigen::internal::set_is_malloc_allowed(true);
  std::cout << "x=" << x.transpose() << std::endl;
  std::cout << "y=" << y.transpose() << std::endl;
  std::cout << "z=" << z.transpose() << std::endl;
}
BOOST_AUTO_TEST_CASE(testAllocationFOO3)
{
  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::VectorXd y = Eigen::VectorXd::Random(4);
  Eigen::VectorXd z(4);
  Eigen::internal::set_is_malloc_allowed(false);
  {
    foo3(z, x, y);
  }
  Eigen::internal::set_is_malloc_allowed(true);
  std::cout << "x=" << x.transpose() << std::endl;
  std::cout << "y=" << y.transpose() << std::endl;
  std::cout << "z=" << z.transpose() << std::endl;
}
BOOST_AUTO_TEST_CASE(testAllocationFOO4)
{
  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::VectorXd y = Eigen::VectorXd::Random(4);
  Eigen::VectorXd z(4);
  Eigen::internal::set_is_malloc_allowed(false);
  {
    foo4(z, x, y);
  }
  Eigen::internal::set_is_malloc_allowed(true);
  std::cout << "x=" << x.transpose() << std::endl;
  std::cout << "y=" << y.transpose() << std::endl;
  std::cout << "z=" << z.transpose() << std::endl;
}
