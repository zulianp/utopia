#ifndef UTOPIA_MULTIPHYSICS_HPP
#define UTOPIA_MULTIPHYSICS_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_ProblemBase.hpp"

namespace utopia {

    // values have to be organized to fit a Newton solver
    // f(x) + b; -> J(x) * s = -(f(x) + b), x = x + s

    template <class FunctionSpace>
    class StationaryLinearFEProblem : public ProblemBase<FunctionSpace> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using Super = utopia::ProblemBase<FunctionSpace>;

        virtual ~StationaryLinearFEProblem() = default;

        void read(Input &in) override {
            Super::read(in);
            in.get("assembly", *assembler_);
            output_path_ = this->output_dir() + ("/" + this->name() + ".e");
            in.get("solver", *linear_solver_);
        }

        StationaryLinearFEProblem(const std::shared_ptr<FunctionSpace> &space)
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

    template <class FunctionSpace>
    class EulerLinearFEProblem : public TimeDependentProblem<FunctionSpace> {
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
            in.get("mass_matrix_solver", *mass_matrix_linear_solver_);
            in.get("implicit", implicit_);
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
                utopia::out() << "Implicit = " << implicit_ << '\n';
                utopia::out() << "--------------------------------\n";
            }
        }

        EulerLinearFEProblem(const std::shared_ptr<FunctionSpace> &space)
            : Super(space),
              assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_assembler_(std::make_shared<OmniAssembler_t>(space)),
              linear_solver_(std::make_shared<OmniLinearSolver_t>()),
              mass_matrix_linear_solver_(std::make_shared<OmniLinearSolver_t>()),
              io_(std::make_shared<IO_t>(*space)) {}

        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            assembler_->set_environment(env);
            mass_matrix_assembler_->set_environment(env);
        }

        bool is_linear() const override { return true; }

        bool init() override {
            if (!Super::init()) return false;

            if (!assembler_->assemble(*this->solution(), *this->jacobian(), *this->fun())) {
                return false;
            }

            *this->fun() *= -1.0;
            forcing_function_ = *this->fun();

            assert(Scalar_t(norm2(*this->solution())) == 0.0);

            // material_matrix_ = *this->jacobian();

            mass_matrix_ = std::make_shared<Matrix_t>();
            this->space()->create_matrix(*mass_matrix_);

            this->space()->create_vector(increment_);

            // FIXME remove increment
            mass_matrix_assembler_->assemble(*this->solution(), *mass_matrix_, increment_);
            rename(this->name() + "_mass_matrix", *mass_matrix_);
            increment_.set(0.0);

            assert(this->delta_time() > 0);

            if (implicit_) {
                assert(this->delta_time() > 0);
                (*this->jacobian()) *= this->delta_time();
                (*this->jacobian()) += *mass_matrix_;

                linear_solver_->update(this->jacobian());
            } else {
                (*this->jacobian()) *= -this->delta_time();
                (*this->jacobian()) += *mass_matrix_;

                this->space()->apply_constraints(*mass_matrix_);
                mass_matrix_linear_solver_->update(mass_matrix_);
            }

            // Apply BC so that we can use the increment with zero BC after
            this->space()->apply_constraints(*this->solution());

            this->export_tensors();
            return true;
        }

        bool update() override { return assemble() && apply_constraints() && solve(); }

        /// Being a linear problem the Jacobian is assembled only once in the init function
        bool assemble() override {
            if (implicit_) {
                (*this->fun()) = -((*this->jacobian()) * (*this->solution()));
                (*this->fun()) += this->delta_time() * forcing_function_;
                (*this->fun()) += (*mass_matrix_) * (*this->solution());
            } else {
                *this->fun() = (*this->jacobian()) * (*this->solution());
            }

            return true;
        }

        bool apply_constraints() override {
            if (implicit_) {
                // FIXME this is not efficient if it is linear
                this->space()->apply_constraints(*this->jacobian());
                this->space()->apply_zero_constraints((*this->fun()));
            } else {
                this->space()->apply_constraints(*this->jacobian(), *this->fun());
            }

            return true;
        }

        bool solve() override {
            if (implicit_) {
                increment_.set(0.0);
                bool ok = linear_solver_->apply(*this->fun(), increment_);
                (*this->solution()) += increment_;

                if (verbose_) {
                    utopia::out() << "Step:\t" << this->current_time().step() << "\tTime:\t"
                                  << this->current_time().get() << "\tSum increment " << Scalar_t(sum(increment_))
                                  << '\n';
                }

                return ok;
            } else {
                return mass_matrix_linear_solver_->apply(*this->fun(), *this->solution());
            }
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

        std::vector<std::shared_ptr<Matrix_t>> operators() override {
            std::vector<std::shared_ptr<Matrix_t>> ret{this->jacobian(), mass_matrix_};
            return ret;
        }

    private:
        std::shared_ptr<OmniAssembler_t> assembler_;
        std::shared_ptr<OmniAssembler_t> mass_matrix_assembler_;
        Path output_path_;
        std::shared_ptr<LinearSolver_t> linear_solver_;
        std::shared_ptr<LinearSolver_t> mass_matrix_linear_solver_;
        std::shared_ptr<Matrix_t> mass_matrix_;
        bool implicit_{false};
        bool verbose_{false};

        std::shared_ptr<IO_t> io_;
        Vector_t increment_, forcing_function_;
        // Matrix_t material_matrix_;
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
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using StationaryLinearFEProblem_t = utopia::StationaryLinearFEProblem<FunctionSpace>;
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

        bool condense() {
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
        using StationaryLinearFEProblem_t = utopia::StationaryLinearFEProblem<FunctionSpace>;
        using EulerLinearFEProblem_t = utopia::EulerLinearFEProblem<FunctionSpace>;
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

                    if (type == "EulerLinearFEProblem") {
                        p = std::make_shared<EulerLinearFEProblem_t>(s);
                    } else {
                        // Default stationary/linear problem
                        p = std::make_shared<StationaryLinearFEProblem_t>(s);
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

        bool reassemble_problems() {
            for (auto &p : problems) {
                bool ok = p.second->assemble();

                if (!ok) {
                    utopia::err() << "Problem between " << p.second->name() << " initialization failed!\n";
                    Utopia::Abort();
                }
            }

            return true;
        }

        bool assemble_problems() {
            for (auto &p : problems) {
                bool ok = p.second->init() && p.second->assemble();

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
                if (true) {
                    // if (problems.size() == 1UL) {
                    do {
                        increment_time();
                        if (!(p->update() && export_results())) {
                            return false;
                        }

                    } while (!p->complete());
                }

                // else {
                //     // Branch is buggy

                //     bool trivial = true;
                //     for (auto &p : problems) {
                //         trivial = trivial && p.second->is_trivial();
                //     }

                //     if (trivial) {
                //         increment_time();

                //         do {
                //             if (!(p->solve() &&             //
                //                   transfer_solutions() &&   //
                //                   reassemble_problems() &&  //
                //                   // apply_constraints() &&     //
                //                   // transfer_residuals() &&    //
                //                   p->apply_constraints() &&  //
                //                   export_results())) {
                //                 return false;
                //             }

                //             increment_time();

                //         } while (!p->complete());
                //     } else {
                //         assert(false && "IMPLEMENT ME");
                //         return false;
                //     }
                // }

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
