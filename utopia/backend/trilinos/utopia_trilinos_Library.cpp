#include "utopia_trilinos_Library.hpp"
#include "utopia_Base.hpp"

#include <Tpetra_Core.hpp>

namespace utopia {
    void TrilinosLibrary::init(int argc, char *argv[]) { Tpetra::initialize(&argc, &argv); }

    int TrilinosLibrary::finalize() {
        Tpetra::finalize();
        return 0;
    }

    std::string TrilinosLibrary::name() const { return "Trilinos"; }

    std::string TrilinosLibrary::version_string() const { return "TODO"; }

}  // namespace utopia
