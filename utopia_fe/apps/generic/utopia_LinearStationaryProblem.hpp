#ifndef UTOPIA_LINEAR_STATIONARY_PROBLEM_HPP
#define UTOPIA_LINEAR_STATIONARY_PROBLEM_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_ProblemBase.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

namespace utopia {

    template <class FunctionSpace>
    class LinearStationaryProblem : public ProblemBase<FunctionSpace> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using Super = utopia::ProblemBase<FunctionSpace>;

        virtual ~LinearStationaryProblem() = default;

        void read(Input &in) override {
            Super::read(in);
            in.get("assembly", *assembler_);
            output_path_ = this->output_dir() + ("/" + this->name() + ".e");
            in.get("solver", *linear_solver_);
        }

        LinearStationaryProblem(const std::shared_ptr<FunctionSpace> &space)
            : Super(space),
              assembler_(std::make_shared<OmniAssembler_t>(space)),
              linear_solver_(std::make_shared<OmniLinearSolver_t>()) {}

        void set_environment(const std::shared_ptr<Environment_t> &env) override { assembler_->set_environment(env); }

        bool is_linear() const override { return true; }

        bool init() override {
            if (!Super::init()) return false;

            if (!assembler_->assemble(*this->solution(), *this->jacobian(), *this->fun())) {
                return false;
            }

            (*this->fun()) *= -1.0;
            return true;
        }

        /// Being a linear problem it is assembled only once in the init function
        bool assemble() override { return true; }

        bool solve() override { return linear_solver_->solve(*this->jacobian(), *this->fun(), *this->solution()); }
        bool export_result() const override { return this->space()->write(output_path_, *this->solution()); }

    private:
        std::shared_ptr<OmniAssembler_t> assembler_;

        Path output_path_;
        std::shared_ptr<LinearSolver_t> linear_solver_;
    };

}  // namespace utopia
#endif  // UTOPIA_LINEAR_STATIONARY_PROBLEM_HPP