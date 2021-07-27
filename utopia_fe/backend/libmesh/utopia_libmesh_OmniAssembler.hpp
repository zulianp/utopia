#ifndef UTOPIA_LIBMESH_OMNI_ASSEMBLER_HPP
#define UTOPIA_LIBMESH_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace_new.hpp"

namespace utopia {
    namespace libmesh {

        class OmniAssembler : public Configurable {
        public:
            using Matrix = Traits<libmesh::FunctionSpace>::Matrix;
            using Vector = Traits<libmesh::FunctionSpace>::Vector;

            OmniAssembler(const std::shared_ptr<libmesh::FunctionSpace> &space);
            virtual ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun);
            bool assemble(const Vector &x, Matrix &jacobian);
            bool assemble(const Vector &x, Vector &fun);
            void read(Input &in) override;

            void set_environment(const std::shared_ptr<Environment<libmesh::FunctionSpace>> &env);

            bool is_linear() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace libmesh

    template <>
    class OmniAssembler<utopia::libmesh::FunctionSpace> final : public utopia::libmesh::OmniAssembler {
    public:
        using Super = utopia::libmesh::OmniAssembler;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_OMNI_ASSEMBLER_HPP
