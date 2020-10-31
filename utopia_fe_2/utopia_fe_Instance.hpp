#ifndef UTOPIA_FE_INSTANCE_HPP
#define UTOPIA_FE_INSTANCE_HPP

#include "utopia_Library.hpp"

#include <memory>
#include <vector>

namespace utopia {
    class UtopiaFE {
    public:
        static void Init(int argc, char *argv[]);
        static int Finalize();

        static UtopiaFE &instance();
        inline void add_library(std::unique_ptr<Library> &&l) { libraries_.push_back(std::move(l)); }

    private:
        std::vector<std::unique_ptr<Library>> libraries_;
    };
}  // namespace utopia
#endif  // UTOPIA_FE_INSTANCE_HPP