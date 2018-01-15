// This demo program solves a hyperelastic problem

#include <dolfin.h>
#include "LeastSquaresElast.h"
#include <utopia.hpp>
#include "../FenicsUtopiaFunction.hpp"
#include "../FenicsUtopiaGLFunction.hpp"

using namespace dolfin;

// Sub domain for clamp at left end
class Left : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (std::abs(x[0]) < DOLFIN_EPS) && on_boundary;
  }
};

// Sub domain for Press at right end
class Right : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (std::abs(x[0] - 1.0) < DOLFIN_EPS) && on_boundary;
  }
};

// top of cube 
class Up : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (std::abs(x[1]) < DOLFIN_EPS) && on_boundary;
  }
};

// bottom of cube 
class Down : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (std::abs(x[1] - 1.0) < DOLFIN_EPS) && on_boundary;
  }
};

// Dirichlet boundary condition for clamp at left end
class Clamp : public Expression
{
public:
  Clamp() : Expression(2) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.2;
    values[1] = 0.0;
  }
};

// Dirichlet boundary condition for Press at right end
class Press : public Expression
{
public:
  Press() : Expression(2) {}
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = -0.2;
    values[1] = -0.1;
  }

};

// Dirichlet boundary condition for Press at right end
class Zero : public Expression
{
public:
  Zero() : Expression(2, 2) {}
  void eval(Array<double>& values, const Array<double>& x) const
  {

  }

};

int main()
{
  // parameters["ghost_mode"] = "shared_facet"; 
  parameters["ghost_mode"] = "none"; 

  // mesh and define function space
  auto mesh = std::make_shared<UnitSquareMesh>(25, 25);
  // auto mesh = std::make_shared<UnitTriangleMesh>();
  auto V    = std::make_shared<LeastSquaresElast::FunctionSpace>(mesh);

  // Dirichlet boundaries
  auto left  = std::make_shared<Left>();
  auto right = std::make_shared<Right>(); 
  auto up    = std::make_shared<Up>();
  auto down  = std::make_shared<Down>();

  // Dirichlet boundary functions
  auto c    = std::make_shared<Clamp>();
  auto r    = std::make_shared<Press>();
  auto zero = std::make_shared<Zero>();


  // boundary condition
  auto bcl = std::make_shared<const dolfin::DirichletBC>(V->sub(0), c, left);
  auto bcr = std::make_shared<const dolfin::DirichletBC>(V->sub(0), r, right);

  auto zero_up = std::make_shared<const dolfin::DirichletBC>(V->sub(1), zero, up);
  auto zero_down = std::make_shared<const dolfin::DirichletBC>(V->sub(1), zero, down);
  std::vector<const dolfin::DirichletBC *> bcs = {{ bcl.get(), bcr.get(), zero_up.get(), zero_down.get() }};

  // source and boundary traction functions
  // auto f = std::make_shared<Constant>(1.0, 1.0);
  auto f = std::make_shared<Constant>(.0, .0);

  // material parameters
  // const double E  = 10.0;
  // const double nu = 0.3;
  // auto mu     = std::make_shared<Constant>(E/(2*(1 + nu)));
  // auto lambda = std::make_shared<Constant>(E*nu/((1 + nu)*(1 - 2*nu)));

  auto mu     = std::make_shared<Constant>(1);
  auto lambda = std::make_shared<Constant>(1);

  // auto a = std::make_shared<LeastSquaresElast::BilinearForm>(V, V);
  auto a = std::make_shared<LeastSquaresElast::BilinearForm>(V, V, mu, lambda);
  auto L = std::make_shared<LeastSquaresElast::LinearForm>(V, f);

  // {
  //   std::vector<std::shared_ptr<const dolfin::DirichletBC> > bcs_ptrs = {{ bcl, bcr }};
  //   dolfin::PetscMatrix A; 
  //   dolfin::PetscVector b;
  //   dolfin::SystemAssembler assembler(a, L, bcs_ptrs);
  //   assembler.assemble(b);
  
  //   assembler.assemble(A);
  //   utopia::DSMatrixd utopia_A;
  //   utopia::convert(A.mat(), utopia_A);
  //   // disp(utopia_A);
  //   // write("A.m", utopia_A);


  //   utopia::DVectord utopia_b;
  //   utopia::convert(b.vec(), utopia_b);
  //   disp(utopia_b);
  // }

  Function u(V);
  solve(*a == *L, u, bcs);

  {
    // utopia::DVectord utopia_u;
    // std::cout << u.vector()->str(true) << std::endl;
    // utopia::convert(u_d, utopia_u);
    // disp(utopia_u);
  }

  // Save solution in VTK format
  File file("LeastSquaresElast.pvd");
  file << u[0];
  file << u[1][0];
  file << u[1][1];
  file << u[1][2];
  file << u[1][3];

  // Plot solution
  // plot(u[0]);
  // interactive();
  return 0;
}
