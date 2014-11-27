#include <iostream>
#include "Point.hpp"
#include "PointOfR.hpp"

using namespace pgs;
int main()
{
  std::cout << "Hello World!";
  Point p = PointOfR(3);
  std::cout << p.dim() << std::endl;
  return 0;
}
