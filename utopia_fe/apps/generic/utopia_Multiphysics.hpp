#ifndef UTOPIA_MULTIPHYSICS_HPP
#define UTOPIA_MULTIPHYSICS_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_MatrixTransformer.hpp"
#include "utopia_ProblemBase.hpp"

#include "utopia_LinearImplicitTimeDependentProblem.hpp"
#include "utopia_LinearStationaryProblem.hpp"

namespace utopia {

    template <class FunctionSpace>
    class LagrangeMultiplier final : public Configurable {
    public:
        void read(Input &in) override { transfer.read(in); }
        FETransfer<FunctionSpace> transfer;
    };

    template <class FunctionSpace>
    class Coupling final : public Configurable {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using LinearStationaryProblem_t = utopia::LinearStationaryProblem<FunctionSpace>;
        using ProblemBase_t = utopia::ProblemBase<FunctionSpace>;

        std::shared_ptr<ProblemBase_t> from;
        std::shared_ptr<ProblemBase_t> to;
        LagrangeMultiplier<FunctionSpace> lagrange_multiplier;

        void read(Input &in) override {
            in.get("verbose", verbose_);
            lagrange_multiplier.read(in);
        }

        void set(const std::shared_ptr<ProblemBase_t> &from, const std::shared_ptr<ProblemBase_t> &to) {
            this->from = from;
            this->to = to;
        }

        bool init() { return lagrange_multiplier.transfer.init(from->space(), to->space()); }

        bool condense_operators() {
            auto &&from_ops = from->operators();
            auto &&to_ops = to->operators();

            const Size_t n = from_ops.size();
            assert(n == to_ops.size());

            if (verbose_) {
                utopia::out() << "Condensing " << n << " pair(s) of systems!\n";
            }

            for (Size_t i = 0; i < n; ++i) {
                Matrix_t temp;
                if (!lagrange_multiplier.transfer.apply(*to_ops[i], temp)) {
                    return false;
                }

                (*from_ops[i]) += temp;
            }

            return true;
        }

        bool condense_right_hand_sides() {
            auto &&from_rhs = from->right_hand_sides();
            auto &&to_rhs = to->right_hand_sides();

            const Size_t n = from_rhs.size();
            assert(n == to_rhs.size());

            if (verbose_) {
                utopia::out() << "Condensing " << n << " pair(s) of rhs!\n";
            }

            for (Size_t i = 0; i < n; ++i) {
                Vector_t temp;
                if (!lagrange_multiplier.transfer.apply_transpose(*to_rhs[i], temp)) {
                    return false;
                }

                (*from_rhs[i]) += temp;
            }

            return true;
        }

        bool condense() {
            if (condense_operators() && condense_right_hand_sides()) {
                from->condensed_system_built();
                return true;
            } else {
                return false;
            }
        }

        bool transfer_solution() { return lagrange_multiplier.transfer.apply(*from->solution(), *to->solution()); }

        bool transfer_residual() {
            Vector_t temp_residual;
            if (!lagrange_multiplier.transfer.apply_transpose(*to->fun(), temp_residual)) {
                return false;
            }

            (*from->fun()) += temp_residual;
            return true;
        }

    private:
        bool verbose_{false};
    };

    template <class FunctionSpace>
    class Multiphysics final : public Configurable {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using LinearStationaryProblem_t = utopia::LinearStationaryProblem<FunctionSpace>;
        using LinearImplicitTimeDependentProblem_t = utopia::LinearImplicitTimeDependentProblem<FunctionSpace>;
        using ProblemBase_t = utopia::ProblemBase<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;

        void read(Input &in) override {
            in.get("spaces", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    auto s = std::make_shared<FunctionSpace>(comm_);

                    bool read_state = false;
                    node.get("read_state", read_state);
                    if (read_state) {
                        auto field = std::make_shared<Field<FunctionSpace>>();
                        s->read_with_state(node, *field);
                        env->add_field(field);

                    } else {
                        s->read(node);
                    }

                    if (s->name().empty()) {
                        utopia::err() << "name must be defined for space node\n";
                        Utopia::Abort();
                    }

                    env->add_space(s);
                });
            });

            in.get("problems", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    std::string type = "";
                    std::string space = "";
                    node.get("type", type);
                    node.get("space", space);

                    if (space.empty()) {
                        utopia::err() << "space field must be defined for problem node\n";
                        Utopia::Abort();
                    }

                    auto s = env->find_space(space);

                    if (!s) {
                        utopia::err() << "No space with name" + space + "\n";
                        Utopia::Abort();
                    }

                    std::shared_ptr<ProblemBase_t> p;

                    if (type == "LinearImplicitTimeDependentProblem") {
                        p = std::make_shared<LinearImplicitTimeDependentProblem_t>(s);
                    } else {
                        // Default stationary/linear problem
                        p = std::make_shared<LinearStationaryProblem_t>(s);
                    }

                    p->set_environment(env);
                    p->read(node);

                    auto ret = problems.insert(std::make_pair(p->name(), p));
                    if (!ret.second) {
                        utopia::err() << "Problem with name " + p->name() +
                                             " already exists, name field must be unique\n";
                        Utopia::Abort();
                    }
                });
            });

            in.get("couplings", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    auto c = utopia::make_unique<Coupling<FunctionSpace>>();
                    c->read(node);

                    std::string from, to;

                    node.get("from", from);
                    node.get("to", to);

                    assert(!from.empty());
                    assert(!to.empty());

                    auto it_from = problems.find(from);
                    if (it_from == problems.end()) {
                        utopia::err() << "No problem with name " + from + " from field must be defined with valid id\n";
                        Utopia::Abort();
                    }

                    auto it_to = problems.find(to);
                    if (it_to == problems.end()) {
                        utopia::err() << "No problem with name " + to + " to field must be defined with valid id\n";
                        Utopia::Abort();
                    }

                    c->set(it_from->second, it_to->second);
                    couplings.push_back(std::move(c));
                });
            });

            in.get("master_problem", master_problem);

            if (master_problem.empty()) {
                utopia::err() << "master_problem field must be defined with valid problem id\n";
                Utopia::Abort();
            }

            auto it_from = problems.find(master_problem);
            if (it_from == problems.end()) {
                utopia::err() << "No problem with name " + master_problem +
                                     ", master_problem must be defined with valid problem id\n";
                Utopia::Abort();
            }
        }

        bool assemble_couplings() {
            for (auto &c : couplings) {
                bool ok = c->init();

                if (!ok) {
                    utopia::err() << "Coupling between " << c->from->name() << " and " << c->to->name() << " failed!\n";
                    Utopia::Abort();
                }
            }

            return true;
        }

        bool reassemble_operators() {
            for (auto &p : problems) {
                bool ok = p.second->assemble_operators();

                if (!ok) {
                    Utopia::Abort("Problem \"" + p.second->name() + "\": initialization failed!\n");
                }
            }

            return true;
        }

        bool assemble_operators() {
            for (auto &p : problems) {
                bool ok = p.second->init() && p.second->assemble_operators();

                if (!ok) {
                    Utopia::Abort("Problem \"" + p.second->name() + "\": initialization failed!\n");
                }
            }

            return true;
        }

        bool prepare_systems() {
            for (auto &p : problems) {
                bool ok = p.second->prepare_system();

                if (!ok) {
                    Utopia::Abort("System preparation for \"" + p.second->name() + "\" failed!\n");
                }
            }

            return true;
        }

        bool condense_problems() {
            for (auto &c : couplings) {
                if (!c->condense()) {
                    Utopia::Abort("Condesation of " + c->from->name() + " and " + c->to->name() + " failed!\n");
                }
            }

            return true;
        }

        bool apply_constraints() {
            auto it = problems.find(master_problem);
            if (it == problems.end()) {
                assert(false && "SHOULD NEVER COME HERE");
                return false;
            }

            it->second->apply_constraints();
            return true;
        }

        bool assemble_all() {
            // FIXME add stage where residual and system is created
            return assemble_couplings() &&  //
                   assemble_operators() &&  //
                   condense_problems() &&   //
                   prepare_systems() &&     //
                   apply_constraints();     //
        }

        void increment_time() {
            for (auto &p : problems) {
                p.second->increment_time();
            }
        }

        bool solve() {
            auto it = problems.find(master_problem);
            if (it == problems.end()) {
                assert(false && "SHOULD NEVER COME HERE");
                return false;
            }

            auto p = it->second;
            if (p->is_time_dependent()) {
                do {
                    increment_time();
                    if (!(p->update() && export_results())) {
                        return false;
                    }

                } while (!p->complete());

                return true;
            } else {
                bool ok = p->solve();
                return export_results() && ok;
            }
        }

        bool transfer_solutions() {
            for (auto it = couplings.rbegin(); it != couplings.rend(); ++it) {
                auto &c = *it;

                if (!c->transfer_solution()) {
                    utopia::err() << "Projection from " << c->from->name() << " to " << c->to->name() << " failed!\n";
                    Utopia::Abort();
                }
            }
            return true;
        }

        bool transfer_residuals() {
            for (auto it = couplings.rbegin(); it != couplings.rend(); ++it) {
                auto &c = *it;

                if (!c->transfer_residual()) {
                    utopia::err() << "Adjoint projection op from " << c->from->name() << " to " << c->to->name()
                                  << " failed!\n";
                    Utopia::Abort();
                }
            }
            return true;
        }

        bool export_results() {
            if (!transfer_solutions()) {
                return false;
            }

            for (auto &p : problems) {
                p.second->export_result();
            }

            return true;
        }

        bool run() {
            if (assemble_all() && solve()) {
                utopia::out() << "[Status] run succesfull!\n";
                return true;
            } else {
                utopia::err() << "[Error] run unsuccesfull!\n";
                return false;
            }
        }

        Multiphysics(const Communicator_t &comm = Communicator_t())
            : comm_(comm), env(std::make_shared<Environment_t>()) {}

    private:
        Communicator_t comm_;
        std::map<std::string, std::shared_ptr<ProblemBase_t>> problems;
        std::vector<std::unique_ptr<Coupling<FunctionSpace>>> couplings;
        std::string master_problem;

        std::shared_ptr<Environment_t> env;
    };
}  // namespace utopia

#endif  // UTOPIA_MULTIPHYSICS_HPP