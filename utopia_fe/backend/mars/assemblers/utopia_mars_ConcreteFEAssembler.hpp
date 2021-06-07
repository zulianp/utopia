#ifndef UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP
#define UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP

#include "utopia_mars_ForwardDeclarations.hpp"

#include "utopia_mars_FEAssembler.hpp"
#include "utopia_mars_FEHandler.hpp"

#include "utopia_mars_FE.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename...>
        class ConcreteFEAssembler : public FEAssembler {
        public:
            using Scalar = typename Traits<FunctionSpace>::Scalar;

            using FEHandler = utopia::mars::FEHandler<DMesh, 1>;
            using DofHandler = typename FEHandler::DofHandler;
            using FEDofMap = typename FEHandler::FEDofMap;
            using SPattern = typename FEHandler::SPattern;

            using FE = utopia::mars::FE<Scalar, SPattern, FEDofMap>;

            virtual ~ConcreteFEAssembler() = default;

            inline void set_space(const std::shared_ptr<FunctionSpace> &space) { space_ = space; }
            inline std::shared_ptr<FunctionSpace> space() const { return space_; }
            // inline std::shared_ptr<MarsMeshType> mars_mesh_ptr() { return space_->mesh().raw_type<MarsMeshType>(); }
            inline std::shared_ptr<FEHandler> handler() { return space_->raw_type<FEHandler>(); }

        private:
            std::shared_ptr<FunctionSpace> space_;
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP
