#ifndef UTOPIA_COUPLED_FE_FUNCTION_HPP
#define UTOPIA_COUPLED_FE_FUNCTION_HPP

#include "utopia_FEFunctionFactory.hpp"
#include "utopia_FEModelFunction.hpp"
#include "utopia_MatrixTransformer.hpp"

namespace utopia {

    template <class FunctionSpace>
    class CoupledFEFunction : public FEFunctionInterface<FunctionSpace> {
    public:
        using Super = utopia::FEFunctionInterface<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using Size_t = typename Traits<FunctionSpace>::SizeType;

        inline static constexpr bool default_verbose() { return true; }

        class FEProblem : public Configurable {
        public:
            std::string name_{"no_name"};
            // std::shared_ptr<FunctionSpace> space_;
            std::shared_ptr<FEFunctionInterface_t> function_;
            std::shared_ptr<Environment<FunctionSpace>> env_;

            bool is_master_{false};

            // Buffers
            std::shared_ptr<Matrix_t> hessian_;
            std::shared_ptr<Matrix_t> mass_matrix_;
            std::shared_ptr<Vector_t> gradient_;
            std::shared_ptr<Vector_t> solution_;

            bool report_solution() { return function_->report_solution(*solution()); }

            inline bool is_linear() const {
                if (function_) {
                    return function_->is_linear();
                }

                return false;
            }

            inline bool is_time_dependent() const {
                if (function_) {
                    return function_->is_time_dependent();
                }

                return false;
            }

            inline bool is_master() const { return is_master_; }
            inline const std::string &name() const { return name_; }

            inline const std::shared_ptr<Matrix_t> &mass_matrix() const {
                assert(mass_matrix_);
                return mass_matrix_;
            }

            inline std::vector<std::shared_ptr<Matrix_t>> matrices() { return {hessian_}; }
            inline std::vector<std::shared_ptr<Vector_t>> vectors() { return {gradient_}; }

            inline std::shared_ptr<Vector_t> solution() {
                ensure_solution();
                return solution_;
            }

            inline std::shared_ptr<Vector_t> gradient() { return gradient_; }
            inline std::shared_ptr<Matrix_t> hessian() { return hessian_; }

            void ensure_hessian_and_gradient() {
                ensure_solution();
                ensure_gradient();
                ensure_hessian();
            }

            void ensure_gradient() {
                if (!gradient_) {
                    gradient_ = std::make_shared<Vector_t>();
                    space()->create_vector(*gradient_);
                }
            }

            void ensure_hessian() {
                if (!hessian_) {
                    hessian_ = std::make_shared<Matrix_t>();
                    space()->create_matrix(*hessian_);
                }
            }

            void ensure_mass_matrix() {
                if (!mass_matrix_) {
                    mass_matrix_ = std::make_shared<Matrix_t>();
                    space()->create_matrix(*mass_matrix_);
                }
            }

            void ensure_solution() {
                if (!solution_) {
                    solution_ = std::make_shared<Vector_t>();
                    space()->create_vector(*solution_);
                    solution_->set(0.0);
                }
            }

            bool assemble_hessian_and_gradient() {
                ensure_hessian_and_gradient();

                if (!function_) {
                    return false;
                }

                return function_->hessian_and_gradient(*solution_, *hessian_, *gradient_);
            }

            bool assemble_hessian() {
                ensure_hessian();
                ensure_solution();

                if (!function_) {
                    return false;
                }

                return function_->hessian(*solution_, *hessian_);
            }

            bool assemble_mass_matrix() {
                ensure_mass_matrix();
                ensure_solution();

                if (!function_) {
                    return false;
                }

                return function_->assemble_mass_matrix(*mass_matrix_);
            }

            bool assemble_gradient() {
                ensure_gradient();
                ensure_solution();

                if (!function_) {
                    return false;
                }

                return function_->gradient(*solution_, *gradient_);
            }

            void assemble_scalar() {}

            void set_environment(const std::shared_ptr<Environment<FunctionSpace>> &env) {
                assert(function_);

                env_ = env;

                if (function_) {
                    function_->set_environment(env);
                }
            }

            void set_function(const std::shared_ptr<FEFunctionInterface_t> &function) {
                function_ = function;
                function_->must_apply_constraints_to_assembled(false);
            }

            std::shared_ptr<FEFunctionInterface<FunctionSpace>> function() const {
                assert(function_);
                return function_;
            }

            const std::shared_ptr<FunctionSpace> &space() const {
                assert(function_->space());
                return function_->space();
            }

            void read(Input &in) override {
                in.require("name", name_);
                in.get("is_master", is_master_);

                assert(env_);

                if (!env_) {
                    Utopia::Abort("CoupledFEFunction::FEProblem::read: set_environment has to be called before read");
                }

                std::shared_ptr<FunctionSpace> space;
                if (!function_ || !this->space()) {
                    std::string space_id;
                    in.require("space", space_id);
                    space = env_->find_space(space_id);

                    if (!space) {
                        Utopia::Abort("CoupledFEFunction::FEProblem::read: space is undefined!");
                    }
                }

                if (!function_) {
                    function_ = FEFunctionFactory<FunctionSpace>::make(space, in);
                    function_->set_environment(env_);
                    function_->must_apply_constraints_to_assembled(false);
                }

                if (function_) {
                    function_->read(in);
                }
            }
        };

        class Coupling : public Configurable {
        public:
            using TransferAssembler = utopia::FETransfer<FunctionSpace>;

            void read(Input &in) override {
                if (transfer_) {
                    transfer_->verbose(verbose_);
                    transfer_->read(in);
                    in.get("verbose", verbose_);
                }
            }

            bool update() { return transfer_->init(from_->space(), to_->space()); }

            void set(const std::shared_ptr<FEProblem> &from, const std::shared_ptr<FEProblem> &to) {
                this->from_ = from;
                this->to_ = to;
            }

            bool condense_mass_matrix() {
                auto &&from_mass_matrix = from_->mass_matrix();
                auto &&to_mass_matrix = to_->mass_matrix();

                if (verbose_) {
                    utopia::out() << "Condensing mass matrix!\n";
                }

                Matrix_t temp;
                if (!transfer_->apply(*to_mass_matrix, temp)) {
                    return false;
                }

                (*from_mass_matrix) += temp;
                return true;
            }

            bool condense_matrices() {
                auto &&from_matrices = from_->matrices();
                auto &&to_matrices = to_->matrices();

                const Size_t n = from_matrices.size();
                assert(n == to_matrices.size());

                if (verbose_) {
                    utopia::out() << "Condensing " << n << " pair(s) of systems!\n";
                }

                for (Size_t i = 0; i < n; ++i) {
                    Matrix_t temp;
                    if (!transfer_->apply(*to_matrices[i], temp)) {
                        return false;
                    }

                    (*from_matrices[i]) += temp;
                }

                return true;
            }

            bool condense_vectors() {
                auto &&from_vectors = from_->vectors();
                auto &&to_vectors = to_->vectors();

                const Size_t n = from_vectors.size();
                assert(n == to_vectors.size());

                if (verbose_) {
                    utopia::out() << "Condensing " << n << " pair(s) of rhs!\n";
                }

                for (Size_t i = 0; i < n; ++i) {
                    Vector_t temp;
                    if (!transfer_->apply_transpose(*to_vectors[i], temp)) {
                        return false;
                    }

                    (*from_vectors[i]) += temp;
                }

                return true;
            }

            bool transfer_solution() { return transfer_->apply(*from_->solution(), *to_->solution()); }

            bool condense_system() {
                if (condense_matrices() && condense_vectors()) {
                    return true;
                } else {
                    return false;
                }
            }

            Coupling() : transfer_(utopia::make_unique<TransferAssembler>()) {}

            inline const std::shared_ptr<FEProblem> &from() const { return from_; }
            inline const std::shared_ptr<FEProblem> &to() const { return to_; }

        private:
            std::unique_ptr<TransferAssembler> transfer_;
            std::shared_ptr<FEProblem> from_, to_;
            bool verbose_{default_verbose()};
        };

        CoupledFEFunction() {}

        virtual ~CoupledFEFunction() = default;

        inline const std::shared_ptr<Matrix_t> &mass_matrix() const override {
            return master_fe_problem_->mass_matrix();
        }

        bool assemble_mass_matrix() override {
            for (auto &fe_pair : fe_problems_) {
                if (!fe_pair.second->assemble_mass_matrix()) {
                    return false;
                }
            }

            for (auto &c : couplings_) {
                if (!c->condense_mass_matrix()) {
                    assert(false);
                    Utopia::Abort("Failed to condense mass matrices!");
                }
            }

            return true;
        }

        bool assemble_mass_matrix(Matrix_t &mass_matrix) override {
            if (assemble_mass_matrix()) {
                mass_matrix = *this->mass_matrix();
                return true;
            } else {
                return false;
            }
        }

        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            for (auto &fe_pair : fe_problems_) {
                fe_pair.second->set_environment(env);
            }
        }

        inline void create_solution_vector(Vector_t &x) override {
            assert(master_fe_problem_);
            master_fe_problem_->function()->create_solution_vector(x);
        }

        inline void apply_constraints(Vector_t &x) const override {
            assert(master_fe_problem_);
            master_fe_problem_->function()->apply_constraints(x);
        }

        bool update(const Vector_t &) override { return true; }

        bool value(const Vector_t &, Scalar_t &v) const override {
            v = -1;
            return false;
        }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            // TODO project x to subproblems
            assert(is_linear());

            (*master_fe_problem_->solution()) = x;
            project_solutions();

            for (auto &fe_pair : fe_problems_) {
                if (!fe_pair.second->assemble_gradient()) {
                    return false;
                }
            }

            for (auto &c : couplings_) {
                if (!c->condense_vectors()) {
                    assert(false);
                    Utopia::Abort("Failed to condense gradients!");
                }
            }

            g = *master_fe_problem_->gradient();

            // master_fe_problem_->function()->gradient(x, g);

            if (must_apply_constraints_) {
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            // TODO project x to subproblems
            assert(is_linear());

            (*master_fe_problem_->solution()) = x;
            project_solutions();

            for (auto &fe_pair : fe_problems_) {
                if (!fe_pair.second->assemble_hessian()) {
                    return false;
                }
            }

            for (auto &c : couplings_) {
                if (!c->condense_matrices()) {
                    assert(false);
                    Utopia::Abort("Failed to condense hessians!");
                }
            }

            H = *master_fe_problem_->hessian();

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
            }

            return true;
        }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            // TODO project x to subproblems
            assert(is_linear());

            (*master_fe_problem_->solution()) = x;
            project_solutions();

            for (auto &fe_pair : fe_problems_) {
                if (!fe_pair.second->assemble_hessian_and_gradient()) {
                    return false;
                }
            }

            for (auto &c : couplings_) {
                if (!c->condense_matrices()) {
                    assert(false);
                    Utopia::Abort("Failed to condense hessians!");
                }

                // if (!c->condense_vectors()) {
                //     assert(false);
                //     Utopia::Abort("Failed to condense gradients!");
                // }
            }

            H = *master_fe_problem_->hessian();
            g = *master_fe_problem_->gradient();

            apply_transformers(H);

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        const std::shared_ptr<FunctionSpace> &space() const override {
            assert(master_fe_problem_);
            return master_fe_problem_->space();
        }

        void read(Input &in) override {
            Super::read(in);

            in.get("verbose", verbose_);

            in.get("spaces", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    std::string name;
                    node.require("name", name);

                    if (env_->find_space(name)) {
                        // Space is already available do not create again
                        return;
                    }

                    auto s = std::make_shared<FunctionSpace>(comm_);

                    bool read_state = false;
                    node.get("read_state", read_state);
                    if (read_state) {
                        auto field = std::make_shared<Field<FunctionSpace>>();
                        s->read_with_state(node, *field);
                        env_->add_field(field);

                    } else {
                        s->read(node);
                    }

                    if (s->name().empty()) {
                        utopia::err() << "name must be defined for space node\n";
                        Utopia::Abort();
                    }

                    env_->add_space(s);
                });
            });

            if (fe_problems_.empty()) {
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

                        auto s = env_->find_space(space);

                        if (!s) {
                            utopia::err() << "No space with name" + space + "\n";
                            Utopia::Abort();
                        }

                        std::shared_ptr<FEProblem> p = std::make_shared<FEProblem>();
                        p->set_environment(env_);
                        p->read(node);

                        auto ret = fe_problems_.insert(std::make_pair(p->name(), p));
                        if (!ret.second) {
                            utopia::err()
                                << "Problem with name " + p->name() + " already exists, name field must be unique\n";
                            Utopia::Abort();
                        }
                    });
                });
            }

            if (couplings_.empty()) {
                in.get("couplings", [this](Input &array_node) {
                    array_node.get_all([this](Input &node) {
                        auto c = utopia::make_unique<Coupling>();
                        c->read(node);

                        std::string from, to;

                        node.get("from", from);
                        node.get("to", to);

                        assert(!from.empty());
                        assert(!to.empty());

                        this->add_coupling(from, to);
                    });
                });
            }

            in.get("matrix_transformers", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    std::string type;
                    node.get("type", type);
                    auto trafo = MatrixTransformer<Matrix_t>::New(type);

                    if (trafo) {
                        transformers_.push_back(std::move(trafo));
                    } else {
                        Utopia::Abort("[Error] Could not find transformer with type " + type + '\n');
                    }
                });
            });

            search_for_master();
        }

        bool is_linear() const override {
            for (auto &ff : fe_problems_) {
                if (!ff.second->is_linear()) {
                    return false;
                }
            }

            return true;
        }

        bool is_time_dependent() const override {
            for (auto &ff : fe_problems_) {
                if (!ff.second->is_time_dependent()) {
                    return false;
                }
            }

            return true;
        }

        inline void must_apply_constraints_to_assembled(const bool val) override { must_apply_constraints_ = val; }

        void add_function(const std::string &name, const std::shared_ptr<FEFunctionInterface_t> &function) {
            auto problem = std::make_shared<FEProblem>();
            problem->function_ = function;
            problem->name_ = name;

            if (env_) {
                problem->set_environment(env_);
            }

            fe_problems_[name] = problem;
        }

        void add_master_function(const std::string &name, const std::shared_ptr<FEFunctionInterface_t> &function) {
            add_function(name, function);

            auto it = fe_problems_.find(name);

            assert(it != fe_problems_.end());

            master_fe_problem_ = it->second;
        }

        void add_coupling(const std::string &from, const std::string &to) {
            auto c = utopia::make_unique<Coupling>();

            auto it_from = fe_problems_.find(from);
            if (it_from == fe_problems_.end()) {
                utopia::err() << "No problem with name " + from + " from field must be defined with valid id\n";
                Utopia::Abort();
            }

            auto it_to = fe_problems_.find(to);
            if (it_to == fe_problems_.end()) {
                utopia::err() << "No problem with name " + to + " to field must be defined with valid id\n";
                Utopia::Abort();
            }

            c->set(it_from->second, it_to->second);
            couplings_.push_back(std::move(c));
        }

        bool report_solution(const Vector_t &x) override {
            *master_fe_problem_->solution() = x;
            project_solutions();

            bool ok = true;
            for (auto &ff : fe_problems_) {
                ok &= ff.second->report_solution();
            }

            return ok;
        }

        void initialize() {
            for (auto &c : couplings_) {
                c->update();
            }
        }

    protected:
        void apply_transformers(Matrix_t &mat) const {
            for (auto &trafo : transformers_) {
                trafo->apply(mat);
            }
        }

    private:
        Communicator_t comm_;
        std::map<std::string, std::shared_ptr<FEProblem>> fe_problems_;
        std::shared_ptr<FEProblem> master_fe_problem_;
        std::vector<std::unique_ptr<Coupling>> couplings_;
        std::shared_ptr<Environment_t> env_;
        bool must_apply_constraints_{true};
        bool verbose_{default_verbose()};

        std::vector<std::unique_ptr<MatrixTransformer<Matrix_t>>> transformers_;

        void project_solutions() const {
            for (auto it = couplings_.rbegin(); it != couplings_.rend(); ++it) {
                auto &c = *it;

                if (!c->transfer_solution()) {
                    utopia::err() << "Projection from " << c->from()->name() << " to " << c->to()->name()
                                  << " failed!\n";
                    Utopia::Abort();
                }
            }
        }

        void search_for_master() {
            for (auto &ff : fe_problems_) {
                if (ff.second->is_master()) {
                    master_fe_problem_ = ff.second;
                    return;
                }
            }

            // If role is not user defined search automatically
            std::set<std::string> names;
            for (auto &ff : fe_problems_) {
                names.insert(ff.second->name());
            }

            // remove all slave problems
            for (auto &c : couplings_) {
                names.erase(c->to()->name());
            }

            if (names.size() == 1) {
                auto master_name = *names.begin();

                for (auto &ff : fe_problems_) {
                    if (ff.second->name() == master_name) {
                        master_fe_problem_ = ff.second;
                        return;
                    }
                }
            }

            assert(false);
            Utopia::Abort("Master problem must be defined! set \"is_master\" key to \"true\" in problem");
        }
    };
}  // namespace utopia

#endif  // UTOPIA_COUPLED_FE_FUNCTION_HPP
