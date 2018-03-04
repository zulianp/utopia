#include "utopia_Instance.hpp"
#include "utopia_Base.hpp"
#include "utopia_Log.hpp"
#include "utopia_Config.hpp"

#ifdef WITH_PETSC
#include "petscsys.h"
#else
#ifdef WITH_MPI
#include <mpi.h>
#endif //WITH_MPI
#endif //WITH_PETSC

namespace utopia {

  void Utopia::Init(int argc, char *argv[]) {

#ifdef WITH_PETSC
    static char help[] = "initializing utopia environment through petsc";

    PetscInitialize(&argc, &argv, (char *) 0, help);
//        PetscInitializeNoArguments();
#else
#ifdef WITH_MPI
    MPI_Init(&argc, &argv);
#endif //WITH_MPI
#endif //WITH_PETSC
  }


  int Utopia::Finalize() {
#ifdef UTOPIA_LOG_ENABLED
    Logger::instance().save_collected_log();
#endif

#ifdef WITH_PETSC
    PetscFinalize();
#else
#ifdef WITH_MPI
    return MPI_Finalize();
#endif //WITH_MPI
#endif //WITH_PETSC

    return 0;
  }

  Utopia &Utopia::Instance()
  {
   static Utopia instance;
   return instance;
 }

 Utopia::Utopia()
 {
  set("data_path", "../data");
  set("opencl_templates_path", "../backend/opencl/templates");
}

bool Utopia::verbose() const
{
  return get("verbose") == "true";
}
}
