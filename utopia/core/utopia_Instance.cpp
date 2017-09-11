#include "utopia_Instance.hpp"
#include "utopia_Base.hpp"
#include "utopia_Log.hpp"
#include "utopia_Config.hpp"

#ifdef WITH_PETSC

#include "petscsys.h"

#endif

namespace utopia {

  void Utopia::Init(int argc, char *argv[]) {

#ifdef WITH_PETSC
    static char help[] = "initializing utopia environment through petsc";

    PetscInitialize(&argc, &argv, (char *) 0, help);
//        PetscInitializeNoArguments();
#endif
  }


  int Utopia::Finalize() {
#ifdef UTOPIA_LOG_ENABLED
    Log::instance().save_collected_log();
#endif

#ifdef WITH_PETSC
    PetscFinalize();
#endif

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
}
