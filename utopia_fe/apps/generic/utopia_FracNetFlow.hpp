#ifndef UTOPIA_FRAC_NET_FLOW_HPP
#define UTOPIA_FRAC_NET_FLOW_HPP

#include "utopia_CoupledFEFunction.hpp"

namespace utopia {

    template <class FunctionSpace>
    class NLSolve : public Configurable {
    public:
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

        using NewtonBase_t = utopia::NewtonBase<Matrix_t, Vector_t>;
        using Newton_t = utopia::Newton<Matrix_t, Vector_t>;

        inline std::shared_ptr<Environment_t> &env() { return env_; }
        inline const Communicator_t &comm() const { return comm_; }
        inline Communicator_t &comm() { return comm_; }

        void init(const std::shared_ptr<FEFunctionInterface_t> &function) { function_ = function; }

        void read(Input &in) override {
            in.get("solver", *solver_);
            in.get("verbose", verbose_);
        }

        bool solve() {
            Vector_t x;
            function_->create_solution_vector(x);

            if (function_->is_linear() && !function_->is_time_dependent()) {
                // Tivial problem, lets keep it simple

                Matrix_t H;  // Hessian/Jacobian
                Vector_t g;  // Postive gradient / Residual
                bool ok = function_->hessian_and_gradient(x, H, g);
                assert(ok);

                Vector_t c(layout(x), 0.0);  // Correction

                // Solve linear problem
                ok = solver_->linear_solver()->solve(H, g, c);
                assert(ok);

                // Correct the solution with respect to the negative gradient
                x -= c;

                function_->report_solution(x);
                return ok;
            } else {
                function_->setup_IVP(x);

                do {
                    if (!solver_->solve(*function_, x)) {
                        utopia::err() << "NLSolve[Error] Solver failed to solve!\n";
                        return false;
                    }

                    function_->update_IVP(x);
                    function_->report_solution(x);
                } while (!function_->is_IVP_solved());
            }

            return true;
        }

        NLSolve() : env_(std::make_shared<Environment_t>()) { init_defaults(); }

        inline bool verbose() const { return verbose_; }

    private:
        Communicator_t comm_;
        std::shared_ptr<Environment_t> env_;
        std::shared_ptr<FEFunctionInterface_t> function_;
        std::shared_ptr<NewtonBase_t> solver_;
        bool verbose_{true};

        void init_defaults() {
            auto linear_solver = std::make_shared<OmniLinearSolver_t>();
            auto newton = std::make_shared<Newton_t>(linear_solver);
            solver_ = newton;
        }
    };

    template <class FunctionSpace>
    class FracNetFlow : public NLSolve<FunctionSpace> {
    public:
        using Super = utopia::NLSolve<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using TimeDependentFunction_t = utopia::TimeDependentFunction<FunctionSpace>;
        using CoupledFEFunction_t = utopia::CoupledFEFunction<FunctionSpace>;

        std::shared_ptr<FEFunctionInterface_t> create_problem(
            Input &in,
            const std::string &problem_type,
            const std::shared_ptr<FunctionSpace> &porous_matrix,
            const std::vector<std::shared_ptr<FunctionSpace>> &fracture_networks) {
            auto problem = std::make_shared<CoupledFEFunction_t>();

            std::shared_ptr<FEFunctionInterface_t> matrix_problem;

            if (fracture_networks.empty()) {
                this->comm().root_print("FracNetFlow[Warning] fracture_networks is undefined");
            }

            in.get("porous_matrix", [&](Input &node) {
                if (this->verbose()) {
                    this->comm().root_print("FracNetFlow[Status] Reading material for: " + porous_matrix->name());
                }

                auto flow = std::make_shared<FEModelFunction_t>(porous_matrix);
                node.require(problem_type, *flow);
                matrix_problem = flow;
                problem->add_master_function(porous_matrix->name(), flow);
            });

            in.get("fracture_networks", [&](Input &array_node) {
                int idx = 0;
                array_node.get_all([&](Input &node) {
                    auto space = fracture_networks[idx++];
                    auto flow = std::make_shared<FEModelFunction_t>(space);
                    node.require(problem_type, *flow);
                    problem->add_function(space->name(), flow);
                    problem->add_coupling(porous_matrix->name(), space->name());
                });
            });

            problem->initialize();

            if (problem_type == "transport") {
                return std::make_shared<ImplicitEulerIntegrator_t>(problem);
            } else {
                return problem;
            }
        }

        void read(Input &in) override {
            Super::read(in);

            std::shared_ptr<FunctionSpace> porous_matrix;
            std::vector<std::shared_ptr<FunctionSpace>> fracture_networks;

            in.get("porous_matrix", [this, &porous_matrix](Input &node) {
                // Read the function-space of the porous-matrix

                auto s = std::make_shared<FunctionSpace>(this->comm());
                // Use this so everyhting is added to the env automatically when calling read
                // s->set_environment(this->env());
                node.require("space", *s);

                if (s->name().empty()) {
                    utopia::err() << "name must be defined for space node\n";
                    Utopia::Abort();
                }

                this->env()->add_space(s);
                porous_matrix = s;
            });

            in.get("fracture_networks", [this, &fracture_networks](Input &array_node) {
                array_node.get_all([this, &fracture_networks](Input &node) {
                    // Read the function-space of the fracture-network

                    auto s = std::make_shared<FunctionSpace>(this->comm());
                    node.require("space", *s);
                    fracture_networks.push_back(s);

                    if (s->name().empty()) {
                        utopia::err() << "name must be defined for space node\n";
                        Utopia::Abort();
                    }

                    this->env()->add_space(s);
                    fracture_networks.push_back(s);
                });
            });

            std::string problem_type = "flow";
            in.get("problem_type", problem_type);

            auto problem = create_problem(in, problem_type, porous_matrix, fracture_networks);

            this->init(problem);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_FRAC_NET_FLOW_HPP
