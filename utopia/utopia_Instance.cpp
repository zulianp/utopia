#include "utopia_Instance.hpp"
#include "utopia_Base.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Config.hpp"

#ifdef WITH_PETSC
#include "petscsys.h"
#include "utopia_petsc_UTOPIA_KSP_Solver.hpp"
#else
#ifdef WITH_MPI
#include <mpi.h>
#endif //WITH_MPI
#endif //WITH_PETSC

#include <cassert>

namespace utopia {
    
    void Utopia::Init(int argc, char *argv[]) {
        
#ifdef WITH_PETSC
        static char help[] = "initializing utopia environment through petsc";
        
        PetscInitialize(&argc, &argv, (char *) 0, help);
        
        // is this proper place for doing this ???
        KSPRegister("utopia", KSPCreate_UTOPIA);
        
        //        PetscInitializeNoArguments();
#else
#ifdef WITH_MPI
        MPI_Init(&argc, &argv);
#endif //WITH_MPI
#endif //WITH_PETSC
    }
    
    
    int Utopia::Finalize() {
#ifdef UTOPIA_TRACE_ENABLED
        Tracer::instance().save_collected_log();
#endif
        
#ifdef WITH_PETSC
        PetscFinalize();
#else
#ifdef WITH_MPI
        return MPI_Finalize();
#endif //WITH_MPI
#endif //WITH_PETSC
        

        instance().logger().flush();
        instance().maintenance_logger().flush();
        return 0;
    }
    
    Utopia &Utopia::instance()
    {
        static Utopia instance;
        return instance;
    }
    
    Utopia::Utopia()
    {
        set("data_path", "../data");
        set("opencl_templates_path", "../backend/opencl/templates");

        logger_ = std::make_shared<StandardLogger>();
      
        auto temp = std::make_shared<StandardLogger>();
        temp->set_direct_output(false, false, false);
        maintenance_logger_ = temp;
    }
    
    bool Utopia::verbose() const
    {
        return get("verbose") == "true";
    }
}
