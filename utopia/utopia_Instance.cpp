#include "utopia_Instance.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_AuthoredWork.hpp"
#include "utopia_Base.hpp"
#include "utopia_CiteUtopia.hpp"
#include "utopia_Config.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_make_unique.hpp"

#ifdef WITH_TRILINOS
#include "utopia_trilinos_Library.hpp"
#endif  // WITH_TRILINOS

#ifdef WITH_PETSC
#include "utopia_petsc_Library.hpp"
#endif  // WITH_PETSC

#ifdef WITH_MPI
#include <mpi.h>
#endif  // WITH_MPI

#include <cassert>
#include <iostream>

namespace utopia {

    void Utopia::Init(int argc, char *argv[]) {
#ifdef WITH_PETSC
        instance().add_library(utopia::make_unique<PetscLibrary>());
#endif

#ifdef WITH_TRILINOS
        instance().add_library(utopia::make_unique<TrilinosLibrary>());
#endif  // WITH_TRILINOS

#ifdef WITH_MPI
        if (instance().libraries_.empty()) {
            MPI_Init(&argc, &argv);
        }
#endif

        for (const auto &l : instance().libraries_) {
            l->init(argc, argv);
        }

        instance().read_input(argc, argv);

        // if (instance().verbose() && mpi_world_rank() == 0) {
        // std::cout << "Available libs:\n";
        // for (const auto &l : instance().libraries_) {
        //     std::cout << "- " << l->name() << "\n";
        // }
        // }

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

        int ret = 0;

        for (const auto &l : instance().libraries_) {
            ret += l->finalize();
        }

#ifdef WITH_MPI
        if (instance().libraries_.empty()) {
            ret += MPI_Finalize();
        }
#endif  // WITH_MPI

        instance().libraries_.clear();

        if (instance().exit_code_ != EXIT_SUCCESS) {
            std::cerr << "[Warning] exiting with code: " << instance().exit_code_ << std::endl;
            return instance().exit_code_;
        } else {
            return ret;
        }
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
        instance().set("citations", "false");

        for (int i = 1; i < argc; i++) {
            std::string str(argv[i]);

            if (str == "-verbose") {
                instance().set("verbose", "true");
            }

            if (str == "--no_citations") {
                instance().set("citations", "false");
            }

            if (str == "--citations") {
                instance().set("citations", "true");
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
