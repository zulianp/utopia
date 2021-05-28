#ifndef UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP
#define UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP

#include "utopia_mars_FEAssembler.hpp"

namespace utopia {
    namespace mars {

        template <class MarsMeshType, typename...>
        class ConcreteFEAssembler : public FEAssembler {
        public:
            virtual ~ConcreteFEAssembler() = default;

            inline void set_space(const std::shared_ptr<FunctionSpace> &space) { space_ = space; }
            inline std::shared_ptr<FunctionSpace> space() const { return space_; }
            inline std::shared_ptr<MarsMeshType> mars_mesh_ptr() { return space_->mesh().raw_type<MarsMeshType>(); }

        private:
            std::shared_ptr<FunctionSpace> space_;
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP
