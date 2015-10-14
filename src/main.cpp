#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/CartesianPower.h>
#include <manifolds/Point.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>

using namespace mnf;

class RobotManifold : public mnf::CartesianProduct
{
 public:
  RobotManifold() {}
  virtual ~RobotManifold() {}

 private:
  virtual bool isElementary() const { return true; }

  virtual std::shared_ptr<mnf::Manifold> getNewCopy_() const
  {
    std::shared_ptr<RobotManifold> copy(new RobotManifold(*this));
    return copy;
  }

  std::string description(const std::string& prefix, bool firstCall) const
  {
    std::stringstream ss;
    if (firstCall)
    {
      ss << prefix << "/robot+++++++++++++++++++++++++++++" << std::endl;
      ss << this->description(prefix + "+ ", false);
      ss << prefix << "\\++++++++++++++++++++++++++++++++++" << std::endl;
    }
    else
    {
      for (size_t i = 0; i < subManifolds_.size(); ++i)
      {
        if (subManifolds_[i]->isElementary())
        {
          ss << subManifolds_[i]->description(prefix);
        }
        else
        {
          ss << prefix << "/robot+++++++++++++++++++++++++++++" << std::endl;
          ss << subManifolds_[i]->description(prefix + "+ ", false);
          ss << prefix << "\\++++++++++++++++++++++++++++++++++" << std::endl;
        }
      }
    }
    return ss.str();
  }
};

int main()
{
  CartesianProduct merged;
  RealSpace R8(8);
  RealSpace R3(3);
  S2 s2;
  SO3<ExpMapQuaternion> so3;
  RobotManifold robot1;
  robot1.multiply(R8);
  robot1.display();
  RobotManifold robot2;
  CartesianProduct se3{&R3, &so3};
  robot2.multiply(se3);
  robot2.display();
  merged.multiply(robot1);
  merged.display();
  merged.multiply(robot2);
  merged.multiply(s2);
  merged.multiply(R3);
  merged.multiply(s2);
  merged.multiply(R3);
  merged.display();
  return 0;
}
