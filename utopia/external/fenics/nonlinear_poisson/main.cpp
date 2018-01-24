/*
* @Author: alenakopanicakova
* @Date:   2016-05-07
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-05-29
*/

// This demo illustrates how to use of Fenics for solving a nonlinear
// PDE by using UTOPIA solvers, in this case a nonlinear variant of Poisson's equation,
//
//     - div (1 + u^2) grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = x*sin(y)
//
// and boundary conditions given by
//
//     u(x, y)     = 1  for x = 0
//     du/dn(x, y) = 0  otherwise
//
// This is equivalent to solving the variational problem
//
//    F(u) = ((1 + u^2)*grad(u), grad(v)) - (f, v) = 0

#include <utopia.hpp>
#include <dolfin.h>
#include "NonlinearPoisson.h"
#include "../FenicsUtopiaFunction.hpp"
#include "../FenicsUtopiaGLFunction.hpp"


using namespace dolfin;

// Right-hand side
class Source : public Expression
{
public:

  Source() : Expression() {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = x[0]*sin(x[1]);
  }

};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return std::abs(x[0] - 1.0) < DOLFIN_EPS && on_boundary;
  }
};

int main()
{
  // Create mesh and define function space
  auto mesh = std::make_shared<UnitSquareMesh>(16, 16);
  auto V = std::make_shared<NonlinearPoisson::FunctionSpace>(mesh);


  // Define boundary condition
  auto dirichlet_boundary = std::make_shared<DirichletBoundary>();
  auto g  = std::make_shared<Constant>(1.0);
  auto bc = std::make_shared<const dolfin::DirichletBC>(V, g, dirichlet_boundary);

  // Define source and solution functions
  auto f = std::make_shared<Source>();
  auto u = std::make_shared<Function>(V);
  
  // gradient 
  auto F = std::make_shared<NonlinearPoisson::LinearForm>(V);
  F->u = u;
  F->f = f;

  // hessian
  auto J = std::make_shared<NonlinearPoisson::JacobianForm>(V, V);
  J->u = u;

  // energy 
  auto Pi = std::make_shared<NonlinearPoisson::Form_Pi>(mesh);
  Pi->u = u;
  Pi->f = f;

  std::vector<std::shared_ptr<const dolfin::DirichletBC>> _bcs = {bc};

  //  form initial guess 
  PetscVector uu; 
  assemble(uu, *F); 
  utopia::DVectord x_0; 
  convert(uu.vec(), x_0); 

  utopia::FenicsUtopiaFunction<utopia::DSMatrixd, utopia::DVectord>  fun(u, Pi, F, J, _bcs); 
  utopia::solve(fun, x_0); 

  // Save solution in VTK format
  // File file("nonlinear_poisson.pvd");
  // file << *u;

  // Plot solution
  plot(u);
  interactive();

  return 0;
}
