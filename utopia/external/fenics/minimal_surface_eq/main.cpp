/*
* @Author: alenakopanicakova
* @Date:   2016-05-07
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-06-15
*/

// This demo illustrates how to use of Fenics for solving a minimal surface equation
// .................................................................... WORK IN PROGRESS ...................................................................

#include <utopia.hpp>
#include <dolfin.h>
#include "MinimalSurface.h"
#include "MinimalSurfaceLocal.h"
#include "../FenicsUtopiaFunction.hpp"
#include "../FenicsUtopiaGLFunction.hpp"
#include "../FenicsUtopiaAsynchGLFunction.hpp"

using namespace dolfin;

// Right-hand side
class Source : public Expression
{
public:

  Source() : Expression() {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 10 * (x[0] - 1)* (x[0]);  
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

// Sub domain for Dirichlet boundary condition
class WholeBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return on_boundary;
  }
};

// Sub domain for Dirichlet boundary condition
class OneSide : public SubDomain
{  
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      if(on_boundary)
      {
        if(x[0] < 1E-12 or x[1] < 1E-12)
        {
            return true; 
        }
        else
        {
            return false; 
        }
      }
      else
      {
        return false; 
      }
  }

};

int main()
{
  // Create mesh and define function space
  parameters["ghost_mode"] = "none"; 
  auto mesh = std::make_shared<UnitSquareMesh>(5, 5);
  auto V = std::make_shared<MinimalSurface::FunctionSpace>(mesh);


  auto dirichlet_boundary = std::make_shared<OneSide>();
  auto g  = std::make_shared<Constant>(0.0);
  auto bc = std::make_shared<const dolfin::DirichletBC>(V, g, dirichlet_boundary);
  std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs = {bc};


  auto u = std::make_shared<Function>(V);

  auto k1 = std::make_shared<Constant>(0.5);
  auto k2 = std::make_shared<Constant>(1.0);
  auto f  = std::make_shared<Constant>(1.);

  auto sigma = std::make_shared<Source>();


  // energy 
  auto Pi = std::make_shared<MinimalSurface::Form_Pi>(mesh);
  Pi->u = u;
  Pi->f = f; 
  Pi->sigma = sigma; 
  Pi->k1 = k1; 
  Pi->k2 = k2; 
  
  // gradient 
  auto F = std::make_shared<MinimalSurface::ResidualForm>(V);
  F->u = u;
  F->f = f; 
  F->sigma = sigma; 
  F->k1 = k1; 
  F->k2 = k2; 

  // hessian
  auto J = std::make_shared<MinimalSurface::JacobianForm>(V, V);
  J->u = u;
  J->k1 = k1; 
  J->k2 = k2; 


  // form initial guess 
  PetscVector uu; 
  assemble(uu, *F); 
  utopia::DVectord x_0;
  uu = (*u->vector()); 
  uu = 0.0; 
  Vec up = uu.vec(); 
  utopia::convert(up, x_0); 
  std::cout<<"global size: "<< x_0.size().get(0) << "\n"; 



//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// ........................................................... DOMAIN decomposition .........................................................................
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// auto mesh_local = std::make_shared<Mesh>(MPI_COMM_SELF);
// MeshEditor editor;
// editor.open(*mesh_local, 2, 2);  


// editor.init_cells(mesh->num_cells()); 
// ufc::cell ufc_cell;
// std::vector<std::size_t> cell_data(3);
// int i = 0; 
// for (CellIterator cell(*mesh); !cell.end(); ++cell)
// {   
//     cell_data[0] = cell->entities(0)[0];
//     cell_data[1] = cell->entities(0)[1];
//     cell_data[2] = cell->entities(0)[2];
//     editor.add_cell(cell->index(), cell_data);      
//     //std::cout<<"---- : "<< cell->index() << "   \n"; 
//     i++; 
// }


// i = 0; 
// editor.init_vertices(mesh->num_vertices()); 
// for (VertexIterator vertex(*mesh); !vertex.end(); ++vertex)
// {
//     editor.add_vertex(i, (*vertex).point());
//     i++; 
//  }

// editor.close();



//   auto V_loc = std::make_shared<MinimalSurfaceLocal::FunctionSpace>(mesh_local);
//   auto k1_loc = std::make_shared<Constant>(0.5);
//   auto k2_loc = std::make_shared<Constant>(1.0);
//   auto f_loc  = std::make_shared<Constant>(1.);
//   auto sigma_loc = std::make_shared<Source>();

//   // local BC 
//   auto dirichlet_boundary_loc = std::make_shared<WholeBoundary>();
//   auto g_loc  = std::make_shared<Constant>(0.0);
//   auto bc_loc = std::make_shared<const dolfin::DirichletBC>(V_loc, g_loc, dirichlet_boundary_loc);
//   std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs_loc = {bc_loc};
//   auto u_loc = std::make_shared<Function>(V_loc);

//   // local energy 
//   auto Pi_loc = std::make_shared<MinimalSurfaceLocal::Form_Pi_loc>(mesh_local);
//   Pi_loc->u_loc = u_loc; 
//   Pi_loc->f_loc = f_loc; 
//   Pi_loc->sigma_loc = sigma_loc; 
//   Pi_loc->k1_loc = k1_loc; 
//   Pi_loc->k2_loc = k2_loc; 
  
//   // local gradient 
//   auto F_loc = std::make_shared<MinimalSurfaceLocal::ResidualForm>(V_loc);
//   F_loc->u_loc = u_loc; 
//   F_loc->f_loc = f_loc;  
//   F_loc->sigma_loc = sigma_loc; 
//   F_loc->k1_loc = k1_loc; 
//   F_loc->k2_loc = k2_loc; 

//   // local hessian
//   auto J_loc = std::make_shared<MinimalSurfaceLocal::JacobianForm>(V_loc, V_loc);
//   J_loc->u_loc = u_loc; 
//   J_loc->k1_loc = k1_loc; 
//   J_loc->k2_loc = k2_loc; 



//   // local initial guess 
//   dolfin::PetscVector uu_loc; 
//   assemble(uu_loc, *F_loc); 
//   info(uu_loc); 
//   info(*u_loc->vector()); 
//   *u_loc->vector() = 999.0; 



//   // const Mesh& mesh_tr = *(F_loc->mesh());



//   // dolfin::TensorLayout tl(MPI_COMM_SELF,{}, 0, {}, {});
//   // dolfin::Scalar local_value;
//   // local_value.init(tl);
//   // dolfin::assemble(local_value, *Pi_loc);


//   // dolfin::Scalar global_value;
//   // dolfin::assemble(global_value, *Pi_loc);



//   // std::cout<<"local energy:   "<< local_value.get_scalar_value()<< "  \n"; 
//   // std::cout<<"global energy:   "<< global_value.get_scalar_value()<< "  \n"; 



//   // dolfin::TensorLayout tv(MPI_COMM_SELF,{}, 0, {}, {}); 
//   // dolfin::PetscVector uu_gl; 
//   // uu_gl.init(tv); 
//   // assemble(uu_gl, *F_loc); 


//  //return 0; 


//   // std::cout<<"size of example: "<< x_0.size().get(0) << "\n"; 
//   // auto decomposition = std::make_shared<utopia::FenicsDecomposition<utopia::DSMatrixd, utopia::DVectord, utopia::Matrixd, utopia::Vectord>> (); 
//   // utopia::FenicsUtopiaAsynchGLFunction<utopia::DSMatrixd, utopia::DVectord, utopia::Matrixd, utopia::Vectord>  fun( u, Pi, F, J, bcs, decomposition);                                                      
//   // utopia::gl_solve(fun, x_0, "APTS"); 



//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  utopia::Parameters params; 
  params.tol(1e-7); 
  params.solver_type("TRUST_REGION"); 
  params.verbose(true);  // something is wrong with verbose => TODO: check it out 
  
  utopia::FenicsUtopiaFunction<utopia::DSMatrixd, utopia::DVectord>  fun( u, Pi, F, J, bcs); 
  // utopia::solve(fun, x_0, "trust_region"); 
  utopia::solve(fun, x_0, params); 


 //Save solution in VTK format
  // File file("displacement.pvd");
  // file << *u;

 // Plot solution
 // plot(u);
 // interactive();

  return 0;
}
