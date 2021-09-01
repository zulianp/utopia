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
            using Scalar = Traits<libmesh::FunctionSpace>::Scalar;
            using SimulationTime = utopia::SimulationTime<Scalar>;

            OmniAssembler(const std::shared_ptr<libmesh::FunctionSpace> &space);
            virtual ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun);
            bool assemble(const Vector &x, Matrix &jacobian);
            bool assemble(const Vector &x, Vector &fun);
            bool apply(const Vector &x, Vector &hessian_times_x);

            bool assemble(Matrix &jacobian);
            void read(Input &in) override;

            void set_environment(const std::shared_ptr<Environment<libmesh::FunctionSpace>> &env);

            bool is_linear() const;

            void set_time(const std::shared_ptr<SimulationTime> &time);

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
