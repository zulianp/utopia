#ifndef UTOPIA_LIBMESH_OMNI_ASSEMBLER_HPP
#define UTOPIA_LIBMESH_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace_new.hpp"

namespace utopia {
    namespace libmesh {

        class OmniAssembler : public Configurable {
        public:
            using Matrix = Traits<libmesh::FunctionSpace>::Matrix;
            using Vector = Traits<libmesh::FunctionSpace>::Vector;

            OmniAssembler(const std::shared_ptr<libmesh::FunctionSpace> &space);
            ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun);
            void read(Input &in) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace libmesh

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_OMNI_ASSEMBLER_HPP