/*
* @Author: kopanicakova
* @Date:   2018-04-05 17:59:33
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-04-06 15:45:05
*/
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
    
    #ifdef WITH_SLEPC
        SlepcInitialize(&argc,&argv,(char*)0, help); // calls PetscInitialize inside
    #else        
        PetscInitialize(&argc, &argv, (char *) 0, help);
    #endif    //WITH_SLEPC

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
        #ifdef WITH_SLEPC
            SlepcFinalize(); // calls PetscFinalize inside 
        #else        
            PetscFinalize();
        #endif    //WITH_SLEPC
        
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
