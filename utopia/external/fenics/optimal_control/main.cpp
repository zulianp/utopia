#include <dolfin.h>
#include "OptimalControl.h"
#include "../FenicsUtopiaFunction.hpp"

using namespace dolfin;

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return on_boundary; }
};


class IGuess : public Expression 
{ 
  public:
    IGuess() : Expression() {}

    void eval(Array<double>& values, const Array<double>& coord) const 
    {
      const double x = coord[0];
      values[0] = 0.0; 
    }
};

int main()
{
  auto mesh = std::make_shared<UnitSquareMesh>(50, 50);
  auto V = std::make_shared<OptimalControl::FunctionSpace>(mesh);

  
  auto zero = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  auto bc = std::make_shared<DirichletBC>(V, zero, boundary);

  // Define solution function
  auto u = std::make_shared<Function>(V);

    // Set material parameters
  auto nu  = std::make_shared<Constant>(10e-4); 
  auto delta = std::make_shared<Constant>(6.8); 
  auto betta = std::make_shared<Constant>(6.8);
  auto z = std::make_shared<Constant>(1.0); 


  // Create (linear) form defining (nonlinear) variational problem
  auto F = std::make_shared<OptimalControl::Form_F>(V);
  F->nu = nu; 
  F->delta = delta; 
  F->u = u;
  F->betta = betta; 
  F->z = z;

  // Create Jacobian dF = F' (for use in nonlinear solver).
  auto J = std::make_shared<OptimalControl::Form_J>(V, V);
  J->nu = nu; 
  J->delta = delta; 
  J->u = u;
  J->betta = betta; 
  J->z = z;


  auto Pi = std::make_shared<OptimalControl::Form_Pi>(mesh);
  Pi->nu = nu; 
  Pi->delta = delta; 
  Pi->u = u;
  Pi->betta = betta; 
  Pi->z = z;  



  auto ug = std::make_shared<IGuess>();
  *u = *ug;

  std::vector<std::shared_ptr<const DirichletBC>> bcs = {bc}; 

  
  // Fenics way.... 
  // Set up the non-linear solver
  // auto problem = std::make_shared<NonlinearVariationalProblem>(F, u, bcs, J);
  // NonlinearVariationalSolver solver(problem);
  // solver.parameters["nonlinear_solver"] = "snes";
  // solver.parameters("snes_solver")["linear_solver"] = "lu";
  // solver.parameters("snes_solver")["maximum_iterations"] = 20;
  // solver.parameters("snes_solver")["report"] = true;
  // solver.parameters("snes_solver")["error_on_nonconvergence"] = false;
  // std::pair<std::size_t, bool> out;
  // out = solver.solve();

  
  utopia::FenicsUtopiaFunction<utopia::DSMatrixd, utopia::DVectord>  fun( u, Pi, F, J, bcs); 
  

  auto Mass_form = std::make_shared<OptimalControl::Form_M>(V,V);

  auto Mass_petsc_fenics = std::make_shared<PETScMatrix>();
  assemble(*Mass_petsc_fenics, *Mass_form); 
  Mat Mass_petsc = Mass_petsc_fenics->mat(); 
  utopia::DSMatrixd Mass_utopia; 
  
  // replace with wrap 
  utopia::convert(Mass_petsc, Mass_utopia); 

  utopia::DVectord rsum; 
  Mass_utopia = utopia::diag(sum(Mass_utopia, 1)); 


  // form initial guess 
  dolfin::PETScVector uu; 
  assemble(uu, *F); 
  utopia::DVectord x_0;
  uu = (*u->vector()); 
  Vec up = uu.vec(); 
  utopia::convert(up, x_0); 

  // auto lsolver = std::make_shared< utopia::Factorization<utopia::DSMatrixd, utopia::DVectord> >();
  // utopia::Newton<utopia::DSMatrixd, utopia::DVectord> nlsolver(lsolver);
  // nlsolver.verbose(true); 
  // nlsolver.atol(1e-9); 
  // nlsolver.solve(fun, x_0); 

  // auto subproblem = std::make_shared<utopia::SteihaugToint<utopia::DSMatrixd, utopia::DVectord> >();
  // subproblem->atol(1e-11);
  // utopia::TrustRegion<utopia::DSMatrixd, utopia::DVectord> tr_solver(subproblem);
  // tr_solver.verbose(true); 
  // tr_solver.atol(1e-9); 
  // tr_solver.solve(fun, x_0); 


  auto linear_solver = std::make_shared<utopia::Factorization<utopia::DSMatrixd, utopia::DVectord> >();      
  utopia::AffineSimilarity<utopia::DSMatrixd, utopia::DVectord> solver(linear_solver); 


  solver.set_scaling_matrix(utopia::local_identity(local_size(Mass_utopia).get(0), local_size(Mass_utopia).get(1))); 
  solver.set_mass_matrix(Mass_utopia); 
  solver.verbose(true);
  solver.use_m(false); 
  // solver.set_m(-1); 
  solver.atol(1e-9); 
  solver.verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE); 
  solver.solve(fun, x_0); 



  File file("op_control_u.pvd");
  file << *u;

  return 0;
}
