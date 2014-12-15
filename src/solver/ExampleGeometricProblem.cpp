#include <pgsolver/solver/ExampleGeometricProblem.h>

namespace pgs
{
  RealSpace ExampleGeometricProblem::R3 = RealSpace(3);   

  ExampleGeometricProblem::ExampleGeometricProblem()
    : Problem(R3),
      a(1.0),
      b(0.0),
      c(-2.0),
      d(0.0),
      R1(0.5),
      R2(1.0)
  {
    assert( R1<=R2 && "Must have R1<=R2");
  }

  void ExampleGeometricProblem::getTangentLB(RefVec out) const
  {
    assert(out.size() == M().dim() && "wrong size");
    Eigen::Vector3d AbsoluteLB(-10, -10, -10);
    out = AbsoluteLB - M().getConstView<T>(x().value(), 0); // I think we should use a better accessor for the values of x. Maybe an implementation of Point::GetView
  }
  void ExampleGeometricProblem::getTangentUB(RefVec out) const
  {
    assert(out.size() == M().dim() && "wrong size");
    Eigen::Vector3d AbsoluteUB(10, 10, 10);
    out = AbsoluteUB - M().getConstView<T>(x().value(), 0); // I think we should use smth like getview here. 
  }

  void ExampleGeometricProblem::evalObj(double& out) const
  {
    // Minimize the norm of the optim variable
    // x = [x1, x2, x3];
    // f(x) = x1^2+x2^2+x3^2;
    Eigen::Vector3d v = phi_x_z().value();
    out = v.squaredNorm();
  }
  void ExampleGeometricProblem::evalObjGrad(RefMat out) const
  {
    assert(out.rows() == 1 && "wrong rows size");
    assert(out.cols() == M().dim() && "wrong cols size");
    // x = [x1, x2, x3];
    // f(x) = x1^2+x2^2+x3^2;
    // df/dx(x) = [2.x1, 2.x2, 2.x3]^T;

    Eigen::Vector3d v = x().value();
    out << 2*v[0], 2*v[1], 2*v[2]; 
    M().applyDiffMap(out, out, x()[0]);
  }

  void ExampleGeometricProblem::evalLinCstr(RefVec out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.size() == linCstrDim(i) && "wrong size");
    if(i == 0)
    {
      // Constraint 0
      // Point on plan of equation a.x1+b.x2+c.x3+d=0
      Eigen::Vector3d v = phi_x_z().value();
      out << a*v[0] + b*v[1] + c*v[2] + d;
    }
    else if(i==1)
    {
      // Constraint 1
      // This constraint does not have a linear part
    }
  }
  void ExampleGeometricProblem::evalLinCstrGrad(RefMat out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.rows() == linCstrDim(i) && "Wrong rows size");
    assert(out.cols() == M().dim() && "Wrong cols size");

    if(i == 0)
    {
      // Constraint 0
      // Point on plan of equation a.x1+b.x2+c.x3+d=0
      // x = [x1, x2, x3];
      // f(x) = a.x1+b.x2+c.x3+d;
      // df/dx(x) = [a, b, c]^T;
      out << a, b, c; 
      M().applyDiffMap(out, out, x()[0]);
    }
    else if(i==1)
    {
      // Constraint 1
      // This constraint does not have a linear part
    }

  }
  void ExampleGeometricProblem::getLinCstrLB(RefVec out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.size() == linCstrDim(i) && "wrong size");
    if(i == 0)
    {
      out << 0;
    }
    else if(i==1)
    {
    }
  }
  void ExampleGeometricProblem::getLinCstrUB(RefMat out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.size() == linCstrDim(i) && "wrong size");
    if(i == 0)
    {
      out << 0;
    }
    else if(i==1)
    {
    }
  }


  void ExampleGeometricProblem::evalNonLinCstr(RefVec out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.size() == nonLinCstrDim(i) && "wrong size");
    if(i == 0)
    {
      // Constraint 0
      // This constraint does not have a nonlinear part
    }
    else if(i==1)
    {
      // Constraint 1
      // Point between spheres of radius R1 and R2
      // R1 < x1^2 + x2^2 + x3^2 < R2
      Eigen::Vector3d v = phi_x_z().value();
      out << v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    }
  }
  void ExampleGeometricProblem::evalNonLinCstrGrad(RefMat out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.rows() == nonLinCstrDim(i) && "wrong row size");
    assert(out.cols() == M().dim() && "Wrong cols size");
    if(i == 0)
    {
      // Constraint 0
      // This constraint does not have a nonlinear part
    }
    else if(i==1)
    {
      // Constraint 1
      // Point between spheres of radius R1 and R2
      // R1 < x1^2 + x2^2 + x3^2 < R2
      Eigen::Vector3d v = x().value();
      out << 2*v[0], 2*v[1], 2*v[2]; 
      M().applyDiffMap(out, out, x().value());
    }
  }
  void ExampleGeometricProblem::getNonLinCstrLB(RefVec out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.size() == nonLinCstrDim(i) && "wrong size");
    if(i == 0)
    {
      // Constraint 0
      // This constraint does not have a nonlinear part
    }
    else if(i==1)
    {
      out << R1*R1;
    }
  }
  void ExampleGeometricProblem::getNonLinCstrUB(RefVec out, size_t i) const
  {
    assert( i<numberOfCstr() && "This constraint index is out of bounds");
    assert(out.size() == nonLinCstrDim(i) && "wrong size");
    if(i == 0)
    {
      // Constraint 0
      // This constraint does not have a nonlinear part
    }
    else if(i==1)
    {
      out << R2*R2;
    }
  }

  size_t ExampleGeometricProblem::numberOfCstr() const
  {
    size_t nb = 2;
    return nb;
  }

  Index ExampleGeometricProblem::linCstrDim(size_t i) const
  {
    std::vector<Index> linDim(2);
    linDim[0] = 1;
    linDim[1] = 0;
    return linDim[i];
  }

  Index ExampleGeometricProblem::nonLinCstrDim(size_t i) const
  {
    std::vector<Index> nonLinDim(2);
    nonLinDim[0] = 0;
    nonLinDim[1] = 1;
    return nonLinDim[i];
  }
}
