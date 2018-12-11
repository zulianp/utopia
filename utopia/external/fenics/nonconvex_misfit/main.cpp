#include <dolfin.h>
#include "NonConvexMisfit.h"
#include "../FenicsUtopiaFunction.hpp"

using namespace dolfin;

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return on_boundary; }
};

class C_Source_map : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = sin(6.*3.14*x[0])*sin(2.*3.14*x[1]); 
  }
};


class IGuess : public Expression 
{ 
  public:
    IGuess(MPI_Comm comm) : Expression(2) 
    {
       dolfin::seed(2 + dolfin::MPI::rank(comm)); 
    }

    void eval(Array<double>& values, const Array<double>& x) const 
    {
      values[0] = -1.5; 
      values[1] = -0.05; 
    }
};

int main()
{
  auto mesh = std::make_shared<UnitSquareMesh>(1000, 1000);
  auto V = std::make_shared<NonConvexMisfit::FunctionSpace>(mesh);

  
  auto zero = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();

  auto bc1 = std::make_shared<DirichletBC>(V->sub(0), zero, boundary);
  auto bc2 = std::make_shared<DirichletBC>(V->sub(1), zero, boundary);
  std::vector<std::shared_ptr<const DirichletBC>> bcs = { }; 


  // Define solution function
  auto u = std::make_shared<Function>(V);
  auto c_to_map_against = std::make_shared<C_Source_map>();


  auto F = std::make_shared<NonConvexMisfit::Form_F>(V);
  F->u=u;
  F->c_to_map_against = c_to_map_against;   

  
  auto J = std::make_shared<NonConvexMisfit::Form_J>(V, V);
  J->u=u; 

  auto H = std::make_shared<NonConvexMisfit::Form_H>(V, V);
  H->u=u; 

  auto Pi = std::make_shared<NonConvexMisfit::Form_Pi>(mesh);
  Pi->c_to_map_against = c_to_map_against; 
  Pi->u=u;


  auto ug = std::make_shared<IGuess>(mesh->mpi_comm());
  *u = *ug;
  
  
  // Fenics way - seems like all snes solvers diverge, as system is non-convex... 
  // auto problem = std::make_shared<NonlinearVariationalProblem>(F, u, bcs, J);
  // NonlinearVariationalSolver solver(problem);
  // solver.parameters["nonlinear_solver"] = "snes";
  // solver.parameters("snes_solver")["linear_solver"] = "lu";
  // solver.parameters("snes_solver")["maximum_iterations"] = 100;
  // solver.parameters("snes_solver")["report"] = true;
  // solver.parameters("snes_solver")["error_on_nonconvergence"] = false;
  // std::pair<std::size_t, bool> out;
  // out = solver.solve();


  // form initial guess 
  dolfin::PETScVector uu; 
  assemble(uu, *F); 
  utopia::DVectord x_0;
  uu = (*u->vector()); 
  Vec up = uu.vec(); 
  utopia::convert(up, x_0); 

  File file1("non_convex_ig.pvd");
  file1 << (*u)[0];
  file1 << (*u)[1];


  utopia::FenicsUtopiaFunction<utopia::DSMatrixd, utopia::DVectord>  fun( u, Pi, F, J, bcs); 


  auto subproblem = std::make_shared<utopia::Lanczos<utopia::DSMatrixd, utopia::DVectord> >();
  // subproblem->pc_type("lu"); 
  subproblem->atol(1e-11);
  utopia::TrustRegion<utopia::DSMatrixd, utopia::DVectord> tr_solver(subproblem);
  tr_solver.verbose(true); 
  tr_solver.atol(1e-8); 
  tr_solver.rtol(1e-9); 
  tr_solver.stol(1e-12); 
  tr_solver.max_it(2000); 
  tr_solver.solve(fun, x_0); 



//   auto Mass_form = std::make_shared<NonConvexMisfit::Form_M>(V,V);

//   auto Mass_petsc_fenics = std::make_shared<PETScMatrix>();
//   assemble(*Mass_petsc_fenics, *Mass_form); 
//   Mat Mass_petsc = Mass_petsc_fenics->mat(); 
//   utopia::DSMatrixd Mass_utopia; 
  
//   // replace with wrap 
//   utopia::convert(Mass_petsc, Mass_utopia); 

//   utopia::DVectord rsum; 
//   Mass_utopia = utopia::diag(sum(Mass_utopia, 1)); 


//   auto linear_solver = std::make_shared<utopia::Factorization<utopia::DSMatrixd, utopia::DVectord> >();      
//   utopia::AffineSimilarity<utopia::DSMatrixd, utopia::DVectord> solver(linear_solver); 


// //  solver.set_scaling_matrix(utopia::local_identity(local_size(Mass_utopia).get(0), local_size(Mass_utopia).get(1))); 
//   solver.set_mass_matrix(Mass_utopia); 
//   solver.verbose(true);
//   // solver.use_m(true); 
//   // solver.set_m(-10); 
//   solver.atol(1e-9); 
//   solver.max_it(2000); 
//   solver.verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE); 
//   solver.solve(fun, x_0); 


  File file("non_convex_c.pvd");
  file << (*u)[0];

  File file2("non_convex_gamma.pvd");
  file2 << (*u)[1];

  return 0;
}
