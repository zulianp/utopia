#include "utopia_fe_Instance.hpp"

#include "utopia_Instance.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_fe_Base.hpp"

#ifdef WITH_LIBMESH
#include "utopia_libmesh_Library.hpp"
#endif  // WITH_LIBMESH

namespace utopia {

    UtopiaFE &UtopiaFE::instance() {
        static UtopiaFE instance_;
        return instance_;
    }

    void UtopiaFE::Init(int argc, char *argv[]) {
#ifdef WITH_LIBMESH
        instance().add_library(utopia::make_unique<LibMeshLibrary>());
#endif

        // Init algebra first
        Utopia::Init(argc, argv);

        // Init fe libraries next
        for (const auto &l : instance().libraries_) {
            l->init(argc, argv);
        }
    }

    int UtopiaFE::Finalize() {
        for (const auto &l : instance().libraries_) {
            l->finalize();
        }

        instance().libraries_.clear();

        return Utopia::Finalize();
    }
}  // namespace utopia
