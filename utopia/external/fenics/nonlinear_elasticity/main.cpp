// This demo program solves a hyperelastic problem

#include <dolfin.h>
#include "HyperElasticity.h"
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

// Sub domain for rotation at right end
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


// front of cube 
class Front : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (std::abs(x[2]) < DOLFIN_EPS) && on_boundary;
  }
};

// back of cube 
class Back : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (std::abs(x[2] - 1.0) < DOLFIN_EPS) && on_boundary;
  }
};


// Dirichlet boundary condition for clamp at left end
class Clamp : public Expression
{
public:

  Clamp() : Expression(3) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
  }

};

// Dirichlet boundary condition for rotation at right end
class Rotation : public Expression
{
public:

  Rotation() : Expression(3) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    const double scale = 0.5;

    // Center of rotation
    const double y0 = 0.5;
    const double z0 = 0.5;

    // Large angle of rotation (60 degrees)
    double theta = 1.04719755;

    // New coordinates
    double y = y0 + (x[1] - y0)*cos(theta) - (x[2] - z0)*sin(theta);
    double z = z0 + (x[1] - y0)*sin(theta) + (x[2] - z0)*cos(theta);

    // Rotate at right end
    values[0] = 0.0;
    values[1] = scale*(y - x[1]);
    values[2] = scale*(z - x[2]);
  }

};

int main()
{
  // parameters["ghost_mode"] = "shared_facet"; 
  parameters["ghost_mode"] = "none"; 


  // mesh and define function space
  auto mesh = std::make_shared<UnitCubeMesh>(2, 2, 2);
  auto V = std::make_shared<HyperElasticity::FunctionSpace>(mesh);

  // Dirichlet boundaries
  auto left = std::make_shared<Left>();
  auto right = std::make_shared<Right>(); 
  auto front = std::make_shared<Front>();
  auto back = std::make_shared<Back>();
  auto up = std::make_shared<Up>();
  auto down = std::make_shared<Down>();

  // Dirichlet boundary functions
  auto c = std::make_shared<Clamp>();
  auto r = std::make_shared<Rotation>();

  // boundary condition
  auto bcl = std::make_shared<const dolfin::DirichletBC>(V, c, left);
  auto bcr = std::make_shared<const dolfin::DirichletBC>(V, r, right);
  std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs = {{bcl, bcr}};


  // Dirichlet boundary functions
  // auto c = std::make_shared<Clamp>();

  // // // boundary condition
  // auto bcl = std::make_shared<const dolfin::DirichletBC>(V, c, left);
  // auto bcr = std::make_shared<const dolfin::DirichletBC>(V, c, right);
  // auto bcu = std::make_shared<const dolfin::DirichletBC>(V, c, up);
  // auto bcd = std::make_shared<const dolfin::DirichletBC>(V, c, down);
  // auto bcf = std::make_shared<const dolfin::DirichletBC>(V, c, front);
  // auto bcb = std::make_shared<const dolfin::DirichletBC>(V, c, back);
  // // std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs = {{ bcl, bcr, bcu, bcd, bcf, bcb }};
  // std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs = {{ bcl, bcu, bcd }};

  // source and boundary traction functions
  auto B = std::make_shared<Constant>(0.0, -0.5, 0.0);
  auto T = std::make_shared<Constant>(0.1,  0.0, 0.0);

  // solution function
  auto u = std::make_shared<Function>(V);


  // material parameters
  const double E  = 10.0;
  const double nu = 0.3;
  auto mu = std::make_shared<Constant>(E/(2*(1 + nu)));
  auto lambda = std::make_shared<Constant>(E*nu/((1 + nu)*(1 - 2*nu)));

  // gradient 
  auto F = std::make_shared<HyperElasticity::ResidualForm>(V);
  F->mu = mu;
  F->lmbda = lambda; 
  F->u = u;
  F->B = B; 
  F->T = T;

  // hessian 
  auto J = std::make_shared<HyperElasticity::JacobianForm>(V, V);
  J->u = u;
  J->mu = mu; 
  J->lmbda = lambda; 

  //  energy 
  auto Pi = std::make_shared<HyperElasticity::Form_Pi>(mesh);
  Pi->mu = mu;
  Pi->lmbda = lambda; 
  Pi->u = u;
  Pi->B = B; 
  Pi->T = T;

  // form initial guess 
  PetscVector uu; 
  dolfin::assemble(uu, *F); 
  uu = 0.0;
  Vec u_petsc = uu.vec();
  utopia::DVectord x_0; 
  convert(u_petsc, x_0); 


  //  std::cout<<"size of example: "<< x_0.size().get(0) << "\n"; 
  //  utopia::FenicsUtopiaFunction<utopia::DSMatrixd, utopia::DVectord>  fun( u, Pi, F, J, bcs); 
  //  utopia::solve(fun, x_0, "trust_region"); 
  //  utopia::solve(fun, x_0, "newton"); 
  
  auto decomposition = std::make_shared<utopia::FenicsDecomposition<utopia::DSMatrixd, utopia::DVectord, utopia::Matrixd, utopia::Vectord>> (); 
  utopia::FenicsUtopiaGLFunction<utopia::DSMatrixd, utopia::DVectord, utopia::Matrixd, utopia::Vectord>  fun( u, Pi, F, J, bcs, decomposition);                                                      
  utopia::gl_solve(fun, x_0, "APTS"); 



  //Save solution in VTK format
  File file("displacement.pvd");
  file << *u;

  // Plot solution
  // plot(*u);
  // interactive();

  return 0;
}
