#include "utopia_Base.hpp"

// FIXME: this file causes nvcc to fail
#ifndef KOKKOS_ENABLE_CUDA

#include "utopia_RangeDevice.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_BCConditions.hpp"
#include "utopia_BCSetup.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_MLSteadyState.hpp"
#include "utopia_QP_Poisson.hpp"

#ifdef WITH_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_RedundantQPSolver.hpp"
#endif  // WITH_PETSC

#include <chrono>
#include <cmath>
#include <random>

namespace utopia {

// // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // //
static void bratu_2d_rmtr(Input &in) {
  static const int Dim = 2;
  static const int NVars = 1;

  using Comm = utopia::PetscCommunicator;
  using Mesh = utopia::PetscStructuredGrid<Dim>;
  using Elem = utopia::PetscUniformQuad4;
  using FunctionSpace = utopia::FunctionSpace<Mesh, NVars, Elem>;

  Comm world;

  MPITimeStatistics stats(world);
  stats.start();

  FunctionSpace space;
  space.read(in);
  stats.stop_and_collect("space-creation");

  stats.start();

  MLSteadyState<FunctionSpace,
                QPPoisson<FunctionSpace>,
                AllZeroBC<FunctionSpace>,
                AllZeroIG<FunctionSpace> > time_stepper(space);

  time_stepper.read(in);
  time_stepper.run();


  disp("------------------- here I am --------------");

  stats.stop_collect_and_restart("end");

  space.comm().root_print(std::to_string(space.n_dofs()) + " dofs");
  stats.stop_and_collect("output");
  stats.describe(std::cout);
}

UTOPIA_REGISTER_APP(bratu_2d_rmtr);

}  // namespace utopia

#endif
