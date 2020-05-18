#include "utopia_Instance.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_AuthoredWork.hpp"
#include "utopia_Base.hpp"
#include "utopia_CiteUtopia.hpp"
#include "utopia_Config.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Tracer.hpp"

#ifdef WITH_TRILINOS
#include <Tpetra_Core.hpp>
#endif  // WITH_TRILINOS

#ifdef WITH_PETSC
#include "petscsys.h"
#include "utopia_petsc_build_ksp.hpp"
#else
#ifdef WITH_MPI
#include <mpi.h>
#endif  // WITH_MPI
#endif  // WITH_PETSC

#include <cassert>
#include <iostream>

namespace utopia {

    // #define DISABLE_LOGGER

    void Utopia::Init(int argc, char *argv[]) {
#ifdef WITH_PETSC
        static char help[] = "initializing utopia environment through petsc";

        PetscOptionsSetValue(nullptr, "-on_error_abort", nullptr);

#ifdef WITH_SLEPC
        SlepcInitialize(&argc, &argv, (char *)nullptr, help);  // calls PetscInitialize inside
#else
        PetscInitialize(&argc, &argv, (char *)0, help);
#endif  // WITH_SLEPC

        // is this proper place for doing this ???
        KSPRegister("utopia", KSPCreate_UTOPIA);

        //        PetscInitializeNoArguments();
#else
#ifdef WITH_MPI
        MPI_Init(&argc, &argv);
#endif  // WITH_MPI
#endif  // WITH_PETSC

#ifdef WITH_TRILINOS
        Tpetra::initialize(&argc, &argv);
#endif  // WITH_TRILINOS

        instance().read_input(argc, argv);

        CitationsDB::instance().cite(Cite<Utopia2016Git>::bibtex());
    }

    int Utopia::Finalize() {
#ifdef UTOPIA_TRACE_ENABLED
        Tracer::instance().save_collected_log();
#endif

        instance().logger().flush();

        if (mpi_world_rank() == 0) {
            instance().maintenance_logger().flush();
        }

        if (instance().get("citations") == "true") {
            CitationsDB::print();
        }

#ifdef WITH_TRILINOS
        Tpetra::finalize();
#endif  // WITH_TRILINOS

#ifdef WITH_PETSC
#ifdef WITH_SLEPC
        SlepcFinalize();  // calls PetscFinalize inside
#else
        PetscFinalize();
#endif  // WITH_SLEPC

#else
#ifdef WITH_MPI
        return MPI_Finalize();
#endif  // WITH_MPI
#endif  // WITH_PETSC

        if (instance().exit_code_ != EXIT_SUCCESS) {
            std::cerr << "[Warning] exiting with code: " << instance().exit_code_ << std::endl;
        }

        return instance().exit_code_;
    }

    Utopia &Utopia::instance() {
        static Utopia instance;
        return instance;
    }

    Utopia::Utopia()

    {
        set("data_path", "../data");
        set("opencl_templates_path", "../backend/opencl/templates");

#ifdef DISABLE_LOGGER
        logger_ = std::make_shared<NullLogger>();
        auto temp = std::make_shared<StandardLogger>();
        maintenance_logger_ = temp;
#else
        logger_ = std::make_shared<StandardLogger>();
        auto temp = std::make_shared<StandardLogger>();
        temp->set_direct_output(false, false, false);
        maintenance_logger_ = temp;
#endif
    }

    bool Utopia::verbose() const { return get("verbose") == "true"; }

    void Utopia::Abort() {
        int error_code = -1;
#ifdef WITH_MPI
        MPI_Abort(MPI_COMM_WORLD, error_code);
#else
        exit(error_code);
#endif
    }

    void Utopia::read_input(int argc, char *argv[]) {
        instance().set("citations", "true");

        for (int i = 1; i < argc; i++) {
            std::string str(argv[i]);

            if (str == "-verbose") {
                instance().set("verbose", "true");
            }

            if (str == "--no_citations") {
                instance().set("citations", "false");
            }

#ifdef ENABLE_NO_ALLOC_REGIONS

            if (str == "-on_alloc_violation_abort") {
                Allocations::instance().abort_on_violation(true);
            }

            if (str == "-mute-allocation-ctrl") {
                Allocations::instance().verbose(false);
            }

            if (str == "-data_path") {
                if (++i >= argc) break;

                if (mpi_world_rank() == 0) {
                    std::cout << "data_path: " << argv[i] << std::endl;
                }

                instance().set("data_path", argv[i]);
            }

#endif  // ENABLE_NO_ALLOC_REGIONS

#ifdef UTOPIA_TRACE_ENABLED
            if (str == "-intercept") {
                if (i + 1 < argc) {
                    Tracer::instance().interceptor().expr(argv[i + 1]);
                    Tracer::instance().interceptor().interrupt_on_intercept(true);
                    std::cout << "Added intercept: " << argv[i + 1] << std::endl;
                }

                i++;
            }

            if (str == "--full_trace") {
                Tracer::instance().full_trace(true);
            }
#endif  // UTOPIA_TRACE_ENABLED
        }
    }

    void Utopia::read(Input &is) {
        for (auto &s : settings_) {
            is.get(s.first, s.second);
        }
    }

    void Utopia::print_usage(std::ostream &os) const {
        for (const auto &s : settings_) {
            os << s.first << " : " << s.second << "\n";
        }
    }

}  // namespace utopia
