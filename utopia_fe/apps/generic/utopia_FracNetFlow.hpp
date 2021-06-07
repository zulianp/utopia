#ifndef UTOPIA_FRAC_NET_FLOW_HPP
#define UTOPIA_FRAC_NET_FLOW_HPP

#include "utopia_CoupledFEFunction.hpp"
#include "utopia_NLSolve.hpp"
#include "utopia_QuadraticFEFunction.hpp"

namespace utopia {

    template <class FunctionSpace>
    class FracNetFlow : public NLSolve<FunctionSpace> {
    public:
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;

        using Super = utopia::NLSolve<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using TimeDependentFunction_t = utopia::TimeDependentFunction<FunctionSpace>;
        using CoupledFEFunction_t = utopia::CoupledFEFunction<FunctionSpace>;
        using QuadraticFEFunction_t = utopia::QuadraticFEFunction<FunctionSpace>;

        std::shared_ptr<FEFunctionInterface_t> create_problem(
            Input &in,
            const std::string &problem_type,
            const std::shared_ptr<FunctionSpace> &porous_matrix,
            const std::vector<std::shared_ptr<FunctionSpace>> &fracture_networks) {
            auto problem = std::make_shared<CoupledFEFunction_t>();

            std::shared_ptr<FEFunctionInterface_t> matrix_problem;

            if (fracture_networks.empty()) {
                this->warning("fracture_networks is undefined");
            }

            in.get("porous_matrix", [&](Input &node) {
                this->status("Reading material for " + porous_matrix->name());

                auto flow = std::make_shared<FEModelFunction_t>(porous_matrix);
                flow->set_environment(this->environment());

                node.require(problem_type, *flow);

                // if (simplify_problem_ && flow->is_linear()) {
                //     auto qf = std::make_shared<QuadraticFEFunction_t>(flow);
                //     matrix_problem = qf;
                // } else {
                matrix_problem = flow;
                // }

                problem->add_master_function(porous_matrix->name(), matrix_problem);
            });

            in.get("fracture_networks", [&](Input &array_node) {
                int idx = 0;
                array_node.get_all([&](Input &node) {
                    auto space = fracture_networks[idx++];

                    this->status("Reading material for " + space->name());

                    auto flow = std::make_shared<FEModelFunction_t>(space);
                    flow->set_environment(this->environment());
                    node.require(problem_type, *flow);

                    // if (simplify_problem_ && flow->is_linear()) {
                    //     auto qf = std::make_shared<QuadraticFEFunction_t>(flow);
                    //     problem->add_function(space->name(), qf);
                    // } else {
                    problem->add_function(space->name(), flow);
                    // }

                    problem->add_coupling(porous_matrix->name(), space->name());
                });
            });

            this->status("Initializing coupled problem");
            problem->initialize();

            std::string integrator;
            in.get("integrator", integrator);

            if (problem_type == "transport") {
                problem->add_matrix_transformer(utopia::make_unique<StabilizeTransport<Matrix_t>>());
            }

            if (problem_type == "transport" || integrator == "ImplicitEuler") {
                if (simplify_problem_ && problem->is_linear()) {
                    return std::make_shared<ImplicitEulerIntegrator_t>(
                        std::make_shared<QuadraticFEFunction_t>(problem));
                } else {
                    return std::make_shared<ImplicitEulerIntegrator_t>(problem);
                }
            } else {
                return problem;
            }
        }

        std::string name() const override { return "FracNetFlow"; }

        void read(Input &in) override {
            Super::read(in);

            in.get("simplify_problem", simplify_problem_);

            std::shared_ptr<FunctionSpace> porous_matrix;
            std::vector<std::shared_ptr<FunctionSpace>> fracture_networks;

            in.get("porous_matrix", [this, &porous_matrix](Input &node) {
                // Read the function-space of the porous-matrix

                node.get("space", [this, &porous_matrix](Input &space_node) {
                    auto s = std::make_shared<FunctionSpace>(this->comm());
                    // Use this so everyhting is added to the env automatically when calling read
                    // s->set_environment(this->environment());

                    bool read_state = false;
                    space_node.get("read_state", read_state);
                    if (read_state) {
                        auto field = std::make_shared<Field<FunctionSpace>>();
                        s->read_with_state(space_node, *field);
                        this->environment()->add_field(field);
                    } else {
                        s->read(space_node);
                    }

                    if (s->name().empty()) {
                        utopia::err() << "Value for key \"name\" must be defined for space node\n";
                        Utopia::Abort();
                    }

                    this->environment()->add_space(s);
                    porous_matrix = s;
                });
            });

            if (!porous_matrix) {
                Utopia::Abort("Definition of \"space\" undefined for node porous_matrix!");
            }

            in.get("fracture_networks", [this, &fracture_networks](Input &array_node) {
                array_node.get_all([this, &fracture_networks](Input &node) {
                    // Read the function-space of the fracture-network

                    node.get("space", [this, &fracture_networks](Input &space_node) {
                        auto s = std::make_shared<FunctionSpace>(this->comm());
                        // Use this so everyhting is added to the env automatically when calling read
                        // s->set_environment(this->environment());

                        bool read_state = false;
                        space_node.get("read_state", read_state);
                        if (read_state) {
                            auto field = std::make_shared<Field<FunctionSpace>>();
                            s->read_with_state(space_node, *field);
                            this->environment()->add_field(field);
                        } else {
                            s->read(space_node);
                        }

                        if (s->name().empty()) {
                            utopia::err() << "Value for key \"name\" must be defined for space node\n";
                            Utopia::Abort();
                        }

                        this->environment()->add_space(s);
                        fracture_networks.push_back(s);
                    });
                });
            });

            std::string problem_type = "flow";
            in.get("problem_type", problem_type);

            auto problem = create_problem(in, problem_type, porous_matrix, fracture_networks);

            problem->read(in);

            this->init(problem);
        }

        bool simplify_problem_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_FRAC_NET_FLOW_HPP
