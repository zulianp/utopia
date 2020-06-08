#ifndef UTOPIA_TRILINOSLIBRARY_HPP
#define UTOPIA_TRILINOSLIBRARY_HPP

#include "utopia_Library.hpp"

namespace utopia {
    class TrilinosLibrary : public Library {
    public:
        void init(int argc, char *argv[]) override;
        int finalize() override;

        std::string name() const override;
        std::string version_string() const override;
    };
}  // namespace utopia

#endif  // UTOPIA_TRILINOSLIBRARY_HPP
