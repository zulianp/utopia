#ifndef UTOPIA_KOKKOS_OMNI_ASSEMBLER_HPP
#define UTOPIA_KOKKOS_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_SimulationTime.hpp"
#include "utopia_Traits.hpp"

#include "utopia_FEAssembler.hpp"
#include "utopia_kokkos_FEAssembler.hpp"

#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE>
        class OmniAssembler : public utopia::FEAssembler<FunctionSpace> {
        public:
            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using SimulationTime = utopia::SimulationTime<Scalar>;

            using Environment = utopia::Environment<FunctionSpace>;

            using Intrepid2FEAssembler = utopia::kokkos::FEAssembler<FE>;
            using Intrepid2FEAssemblerPtr = std::shared_ptr<Intrepid2FEAssembler>;

            OmniAssembler(const std::shared_ptr<FunctionSpace> &space);
            virtual ~OmniAssembler();

            void set_domain_fe(const std::shared_ptr<FE> &fe);

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun) override;
            bool assemble(const Vector &x, Matrix &jacobian) override;
            bool assemble(const Vector &x, Vector &fun) override;

            // For linear only
            bool assemble(Matrix &jacobian) override;
            bool assemble(Vector &fun) override;

            bool apply(const Vector &x, Vector &hessian_times_x) override;

            void read(Input &in) override;

            AssemblyMode mode() const;
            void set_mode(AssemblyMode mode);

            std::string name() const override;

            void set_environment(const std::shared_ptr<Environment> &env) override;
            std::shared_ptr<Environment> environment() const override;
            void set_space(const std::shared_ptr<FunctionSpace> &space) override;
            std::shared_ptr<FunctionSpace> space() const override;
            bool is_linear() const override;
            bool is_operator() const override;

            void add_domain_assembler(const Intrepid2FEAssemblerPtr &assembler);
            void fail_if_unregistered(const bool val);

            void set_time(const std::shared_ptr<SimulationTime> &time) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_OMNI_ASSEMBLER_HPP
