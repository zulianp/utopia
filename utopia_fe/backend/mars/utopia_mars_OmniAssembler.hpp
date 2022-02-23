#ifndef UTOPIA_MARS_OMNIASSEMBLER_HPP
#define UTOPIA_MARS_OMNIASSEMBLER_HPP

#include "utopia_mars_FEAssembler.hpp"

namespace utopia {
    namespace mars {

        class OmniAssembler : public FEAssembler {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;

            OmniAssembler(const std::shared_ptr<mars::FunctionSpace> &space);
            virtual ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &mat, Vector &vec) override;
            bool assemble(const Vector &x, Matrix &mat) override;
            bool assemble(const Vector &x, Vector &vec) override;

            bool assemble(Matrix &mat) override;
            void read(Input &in) override;

            void set_environment(const std::shared_ptr<Environment> &env) override;

            bool is_linear() const override;

            void set_space(const std::shared_ptr<FunctionSpace> &space) override;

            std::string name() const override;
            bool apply(const Vector &x, Vector &hessian_times_x) override;
            bool is_operator() const override;
            bool assemble(Vector &fun) override;
            std::shared_ptr<Environment> environment() const override;
            std::shared_ptr<FunctionSpace> space() const override;
            void set_time(const std::shared_ptr<SimulationTime> &time) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace mars

    template <>
    class OmniAssembler<utopia::mars::FunctionSpace> final : public utopia::mars::OmniAssembler {
    public:
        using Super = utopia::mars::OmniAssembler;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_MARS_OMNIASSEMBLER_HPP