#ifndef UTOPIA_LINEAR_IMPLICIT_TIME_DEPENDENT_PROBLEM_HPP
#define UTOPIA_LINEAR_IMPLICIT_TIME_DEPENDENT_PROBLEM_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_MatrixTransformer.hpp"
#include "utopia_ProblemBase.hpp"

namespace utopia {

    // values have to be organized to fit a Newton solver
    // f(x) + b; -> J(x) * s = -(f(x) + b), x = x + s

    template <class FunctionSpace>
    class LinearImplicitTimeDependentProblem : public TimeDependentProblem<FunctionSpace> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using IO_t = utopia::IO<FunctionSpace>;
        using Super = utopia::TimeDependentProblem<FunctionSpace>;

        void read(Input &in) override {
            Super::read(in);
            in.get("assembly", *assembler_);
            output_path_ = this->output_dir() + ("/" + this->name() + ".e");
            io_->set_output_path(output_path_);

            in.get("solver", *linear_solver_);
            in.get("verbose", verbose_);

            bool user_defined_mass = false;
            in.get("mass", [this, &user_defined_mass](Input &node) {
                mass_matrix_assembler_->read(node);
                user_defined_mass = true;
            });

            if (!user_defined_mass) {
                auto params = param_list(param("material", param_list(param("type", "Mass"), param("lumped", true))));
                mass_matrix_assembler_->read(params);
            }

            if (verbose_) {
                utopia::out() << "--------------------------------\n";
                this->describe(utopia::out().stream());
                utopia::out() << "--------------------------------\n";
            }

            in.get("matrix_transformers", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    std::string type;
                    node.get("type", type);
                    auto trafo = MatrixTransformer<Matrix_t>::New(type);

                    if (trafo) {
                        transformers.push_back(std::move(trafo));
                    } else {
                        Utopia::Abort("[Error] Could not find transformer with type " + type + '\n');
                    }
                });
            });
        }

        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            assembler_->set_environment(env);
            mass_matrix_assembler_->set_environment(env);
        }

        bool assemble_mass_matrix() {
            // FIXME remove increment
            mass_matrix_assembler_->assemble(*this->solution(), *mass_matrix_, increment_);
            rename(this->name() + "_mass_matrix", *mass_matrix_);
            increment_.set(0.0);

            assert(this->delta_time() > 0);

            // Apply BC so that we can use the increment with zero BC after
            this->space()->apply_constraints(*this->solution());
            this->export_tensors();
            return true;
        }

        bool assemble_material() {
            assert(Scalar_t(norm2(*this->solution())) == 0.0);
            if (!assembler_->assemble(*this->solution(), *this->jacobian(), *this->fun())) {
                return false;
            }

            *this->fun() *= -1.0;
            return true;
        }

        bool init() override {
            if (!Super::init()) return false;
            assert(this->delta_time() > 0);

            this->space()->create_matrix(*mass_matrix_);
            this->space()->create_vector(increment_);

            if (assemble_material() && assemble_mass_matrix()) {
                // Apply BC so that we can use the increment with zero BC after
                this->space()->apply_constraints(*this->solution());
                this->export_tensors();
                return true;
            } else {
                return false;
            }
        }

        bool update() override { return assemble() && apply_constraints() && solve(); }

        bool assemble() override {
            compute_system();
            compute_residual();
            return true;
        }

        void apply_transformers() {
            for (auto &trafo : transformers) {
                trafo->apply(*this->jacobian());
            }
        }

        void condensed_system_built() override { apply_transformers(); }

        void compute_residual() {
            residual_ = (*this->jacobian()) * (*this->solution());
            residual_ = (*this->fun()) - residual_;
            residual_ *= this->delta_time();
            residual_ += (*mass_matrix_) * (*this->solution());
        }

        void compute_system() {
            system_ = this->delta_time() * (*this->jacobian());
            system_ += (*mass_matrix_);
        }

        bool apply_constraints() override {
            this->space()->apply_constraints(system_);
            this->space()->apply_zero_constraints(residual_);
            return true;
        }

        bool solve() override {
            increment_.set(0.0);

            bool ok = linear_solver_->solve(system_, residual_, increment_);
            assert(ok);

            (*this->solution()) += increment_;

            if (verbose_) {
                utopia::out() << "Step:\t" << this->current_time().step() << "\tTime:\t" << this->current_time().get()
                              << "\tSum increment " << Scalar_t(sum(increment_)) << '\n';
            }

            return ok;
        }

        bool export_result() const override {
            if (this->integrate_all_before_output()) {
                if (this->complete()) {
                    return io_->write(*this->solution());
                }
            } else {
                return io_->write(*this->solution(), this->current_time().step(), this->current_time().get());
            }

            return true;
        }

        inline std::vector<std::shared_ptr<Matrix_t>> operators() override {
            std::vector<std::shared_ptr<Matrix_t>> ret{this->jacobian(), mass_matrix()};
            return ret;
        }

        inline const std::shared_ptr<Matrix_t> &mass_matrix() const { return mass_matrix_; }
        inline bool is_linear() const override { return true; }
        inline bool has_transformers() const { return !transformers.empty(); }

        LinearImplicitTimeDependentProblem(const std::shared_ptr<FunctionSpace> &space)
            : Super(space),
              assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_assembler_(std::make_shared<OmniAssembler_t>(space)),
              linear_solver_(std::make_shared<OmniLinearSolver_t>()),
              mass_matrix_(std::make_shared<Matrix_t>()),
              io_(std::make_shared<IO_t>(*space)) {
            // Make sure transfromers are available
            register_transfomers<Matrix_t>();
        }

    private:
        std::shared_ptr<OmniAssembler_t> assembler_;
        std::shared_ptr<OmniAssembler_t> mass_matrix_assembler_;
        std::shared_ptr<LinearSolver_t> linear_solver_;
        std::shared_ptr<Matrix_t> mass_matrix_;
        std::shared_ptr<IO_t> io_;

        Vector_t increment_;
        Vector_t residual_;
        Matrix_t system_;
        std::vector<std::unique_ptr<MatrixTransformer<Matrix_t>>> transformers;

        Path output_path_;
        bool verbose_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_LINEAR_IMPLICIT_TIME_DEPENDENT_PROBLEM_HPP