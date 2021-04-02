#ifndef UTOPIA_MULTIPHYSICS_HPP
#define UTOPIA_MULTIPHYSICS_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

namespace utopia {

    template <class FunctionSpace>
    class FEProblem : public Configurable {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

        virtual ~FEProblem() = default;

        void read(Input &in) override {
            in.get("name", name_);
            in.get("assembly", *assembler_);
            output_path_ = "./" + name_ + ".e";
            in.get("output", output_path_);
            in.get("solver", *linear_solver_);
            in.get("export_tensors", export_tensors_);
        }

        FEProblem(const std::shared_ptr<FunctionSpace> &space)
            : space_(space),
              assembler_(std::make_shared<OmniAssembler_t>(space)),
              linear_solver_(std::make_shared<OmniLinearSolver_t>()) {}

        std::shared_ptr<FunctionSpace> space() const { return space_; }
        const std::string &name() const { return name_; }

        bool init() {
            solution_ = std::make_shared<Vector_t>();
            fun_ = std::make_shared<Vector_t>();
            jacobian_ = std::make_shared<Matrix_t>();

            space_->create_matrix(*jacobian_);
            space_->create_vector(*solution_);
            space_->create_vector(*fun_);

            rename(name() + "_jacobian", *jacobian_);
            rename(name() + "_solution", *solution_);
            rename(name() + "_fun", *fun_);

            if (!assembler_->assemble(*solution_, *jacobian_, *fun_)) {
                return false;
            }

            (*fun_) *= -1.0;
            return true;
        }

        bool apply_constraints() {
            space_->apply_constraints(*jacobian_, *fun_);

            if (export_tensors_) {
                write("load_" + jacobian_->name() + ".m", *jacobian_);
                write("load_" + fun_->name() + ".m", *fun_);
            }

            return true;
        }

        std::shared_ptr<Matrix_t> jacobian() const { return jacobian_; }
        std::shared_ptr<Vector_t> solution() const { return solution_; }

        bool solve() { return linear_solver_->solve(*jacobian_, *fun_, *solution_); }

        bool export_result() const { return space_->write(output_path_, *solution_); }

    private:
        std::string name_{"no_name"};
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<OmniAssembler_t> assembler_;
        Path output_path_;
        std::shared_ptr<Matrix_t> jacobian_;
        std::shared_ptr<Vector_t> solution_, fun_;
        std::shared_ptr<LinearSolver_t> linear_solver_;

        bool export_tensors_{false};
    };

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
        using FEProblem_t = utopia::FEProblem<FunctionSpace>;

        std::shared_ptr<FEProblem_t> from;
        std::shared_ptr<FEProblem_t> to;
        LagrangeMultiplier<FunctionSpace> lagrange_multiplier;

        void read(Input &in) override { lagrange_multiplier.read(in); }

        void set(const std::shared_ptr<FEProblem_t> &from, const std::shared_ptr<FEProblem_t> &to) {
            this->from = from;
            this->to = to;
        }

        bool init() { return lagrange_multiplier.transfer.init(from->space(), to->space()); }

        bool condense() {
            Matrix_t temp;
            if (!lagrange_multiplier.transfer.apply(*to->jacobian(), temp)) {
                return false;
            }

            (*from->jacobian()) += temp;
            return true;
        }

        bool project_solution() { return lagrange_multiplier.transfer.apply(*from->solution(), *to->solution()); }
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
        using FEProblem_t = utopia::FEProblem<FunctionSpace>;

        void read(Input &in) override {
            in.get("spaces", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    auto s = std::make_shared<FunctionSpace>(comm_);
                    s->read(node);

                    std::string name = "";
                    node.get("name", name);

                    if (name.empty()) {
                        utopia::err() << "name must be defined for space node\n";
                        Utopia::Abort();
                    }

                    auto ret = spaces.insert(std::make_pair(name, s));
                    if (!ret.second) {
                        utopia::err() << "Space with name " + name + " already exists, name field must be unique\n";
                        Utopia::Abort();
                    }
                });
            });

            in.get("problems", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    std::string space = "";
                    node.get("space", space);

                    if (space.empty()) {
                        utopia::err() << "space field must be defined for problem node\n";
                        Utopia::Abort();
                    }

                    auto it = spaces.find(space);
                    if (it == spaces.end()) {
                        utopia::err() << "No space with name" + space + "\n";
                        Utopia::Abort();
                    }

                    auto p = std::make_shared<FEProblem_t>(it->second);
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

        bool assemble_problems() {
            for (auto &p : problems) {
                bool ok = p.second->init();

                if (!ok) {
                    utopia::err() << "Problem between " << p.second->name() << " initialization failed!\n";
                    Utopia::Abort();
                }
            }

            return true;
        }

        bool condense_systems() {
            for (auto &c : couplings) {
                if (!c->condense()) {
                    utopia::err() << "Condesation of " << c->from->name() << " and " << c->to->name() << " failed!\n";
                    Utopia::Abort();
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
            return assemble_couplings() && assemble_problems() && condense_systems() && apply_constraints();
        }

        bool solve() {
            auto it = problems.find(master_problem);
            if (it == problems.end()) {
                assert(false && "SHOULD NEVER COME HERE");
                return false;
            }

            return it->second->solve();
        }

        bool export_results() {
            for (auto it = couplings.rbegin(); it != couplings.rend(); ++it) {
                auto &c = *it;

                if (!c->project_solution()) {
                    utopia::err() << "Projection from " << c->from->name() << " to " << c->to->name() << " failed!\n";
                    Utopia::Abort();
                }
            }

            for (auto &p : problems) {
                p.second->export_result();
            }

            return true;
        }

        bool run() { return assemble_all() && solve() && export_results(); }

        Multiphysics(const Communicator_t &comm = Communicator_t()) : comm_(comm) {}

    private:
        Communicator_t comm_;
        std::map<std::string, std::shared_ptr<FunctionSpace>> spaces;
        std::map<std::string, std::shared_ptr<FEProblem_t>> problems;
        std::vector<std::unique_ptr<Coupling<FunctionSpace>>> couplings;
        std::string master_problem;
    };
}  // namespace utopia

#endif  // UTOPIA_MULTIPHYSICS_HPP
