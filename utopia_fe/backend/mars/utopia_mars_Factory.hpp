#ifndef UTOPIA_MARS_FACTORY_HPP
#define UTOPIA_MARS_FACTORY_HPP

#include "utopia_Input.hpp"

#include "utopia_mars_ForwardDeclarations.hpp"

#include <memory>

namespace utopia {

    namespace mars {

        class Factory {
        public:
            virtual ~Factory() = default;
            virtual std::unique_ptr<FEAssembler> new_assembler(Input &in) const = 0;

            // Add more mars types when necessary
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_FACTORY_HPP