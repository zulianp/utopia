#ifndef UTOPIA_FE_MODEL_FUNCTION_HPP
#define UTOPIA_FE_MODEL_FUNCTION_HPP

#include "utopia_Input.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_ui.hpp"

#include "utopia_LS_Strategy.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_MatrixTransformer.hpp"
#include "utopia_ProblemBase.hpp"

#include "utopia_SimulationTime.hpp"

#include <limits>

namespace utopia {

    template <class FunctionSpace>
    class FEFunctionInterface
        : public Function<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector>,
          public Operator<typename Traits<FunctionSpace>::Vector> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using SimulationTime = utopia::SimulationTime<Scalar_t>;

        virtual ~FEFunctionInterface() = default;

        bool apply(const Vector_t & /*x*/, Vector_t & /*hessian_times_x*/) const override {
            assert(false);
            Utopia::Abort("Implement me!");
            return false;
        }

        Size size() const override {
            assert(false);
            Utopia::Abort("Implement me!");
            return {};
        }

        Size local_size() const override {
            assert(false);
            Utopia::Abort("Implement me!");
            return {};
        }

        Communicator_t &comm() override = 0;
        const Communicator_t &comm() const override = 0;

        virtual bool is_operator() const { return false; }

        virtual void create_solution_vector(Vector_t &x) = 0;
        virtual void apply_constraints(Vector_t &x) const = 0;
        virtual void set_environment(const std::shared_ptr<Environment_t> &env) = 0;
        virtual void post_solve(Vector_t &) {}

        virtual const std::shared_ptr<Matrix_t> &mass_matrix() const = 0;
        virtual bool assemble_mass_matrix() = 0;
        virtual bool assemble_mass_matrix(Matrix_t &mass_matrix) = 0;
        virtual void set_mass_matrix(const std::shared_ptr<Matrix_t> &) {
            assert(false);
            Utopia::Abort("Implement me!");
        }

        virtual const std::shared_ptr<FunctionSpace> &space() const = 0;

        bool initialize_hessian(Matrix_t &H, Matrix_t &) const override {
            space()->create_matrix(H);
            return true;
        }

        virtual bool is_time_dependent() const = 0;
        virtual bool is_linear() const = 0;

        virtual void must_apply_constraints_to_assembled(const bool) {}
        virtual bool report_solution(const Vector_t &) { return true; }

        virtual void initial_guess_for_solver(Vector_t &) {}

        virtual bool update_BVP() { return false; }
        virtual bool update_IVP(const Vector_t &) { return false; }
        virtual bool setup_IVP(Vector_t &) { return false; }
        virtual bool setup_IVP(IO<FunctionSpace> &) {
            Utopia::Abort("setup_IVP(IO<FunctionSpace> &), not implemented!");
            return false;
        }

        virtual bool register_output(IO<FunctionSpace> &) {
            Utopia::Abort("register_output(IO<FunctionSpace> &), not implemented!");
            return false;
        }

        virtual bool update_output(IO<FunctionSpace> &) {
            Utopia::Abort("update_output(IO<FunctionSpace> &), not implemented!");
            return false;
        }

        virtual void set_output_db(const std::shared_ptr<IO<FunctionSpace>> &) {
            Utopia::Abort("set_output_db(IO<FunctionSpace> &), not implemented!");
        }

        virtual bool is_IVP_solved() { return true; }

        virtual bool set_initial_condition(const Vector_t &) {
            assert(false);
            return false;
        }

        virtual bool time_derivative(const Vector_t &x, Vector_t &dfdt) const {
            UTOPIA_UNUSED(x);
            UTOPIA_UNUSED(dfdt);
            return false;
        }

        virtual void set_time(const std::shared_ptr<SimulationTime> &time) = 0;
        virtual std::shared_ptr<SimulationTime> time() const { return nullptr; }

        virtual std::shared_ptr<LSStrategy<Vector_t>> line_search() { return nullptr; }
    };

    template <class FunctionSpace>
    class FEModelFunction final : public FEFunctionInterface<FunctionSpace> {
    public:
        using Super = utopia::FEFunctionInterface<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using SimulationTime = utopia::SimulationTime<Scalar_t>;

        FEModelFunction(const std::shared_ptr<FunctionSpace> &space)
            : space_(space),
              assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_(std::make_shared<Matrix_t>()) {}

        virtual ~FEModelFunction() = default;

        bool is_linear() const override { return assembler()->is_linear(); }
        bool is_hessian_constant() const override { return is_linear(); }

        inline Communicator_t &comm() override { return space_->mesh().comm(); }
        inline const Communicator_t &comm() const override { return space_->mesh().comm(); }

        bool is_operator() const override { return assembler()->is_operator(); }

        Size size() const override { return {space_->n_dofs(), space_->n_dofs()}; }

        Size local_size() const override { return {space_->n_local_dofs(), space_->n_local_dofs()}; }

        void set_time(const std::shared_ptr<SimulationTime> &time) override {
            assert(assembler_);
            if (assembler_) {
                assembler_->set_time(time);
            }
        }

        void init_mass_matrix_assembler() {
            auto params = param_list(param("material",
                                           param_list(param("type", "Mass"),                         //
                                                      param("lumped", true),                         //
                                                      param("n_components", this->space()->n_var())  //
                                                      )));

            this->mass_matrix_assembler()->read(params);
        }

        void read(Input &in) override {
            Super::read(in);
            in.get("assembly", *assembler_);
            in.get("verbose", verbose_);
            in.get("export_tensors", export_tensors_);

            bool user_defined_mass = false;
            in.get("mass", [this, &user_defined_mass](Input &node) {
                this->mass_matrix_assembler()->read(node);
                user_defined_mass = true;
            });

            if (!user_defined_mass) {
                init_mass_matrix_assembler();
            }

            auto space_name = space()->name();

            if (!space_name.empty()) {
                output_path_ = space_name + ".e";
            }

            in.get("output_path", output_path_);
        }

        bool value(const Vector_t & /*point*/, Scalar_t &value) const override {
            // assert(false && "IMPLEMENT ME");
            // value = std::numeric_limits<Scalar_t>::signaling_NaN();
            value = -12345678;
            return false;
        }

        inline const std::shared_ptr<OmniAssembler_t> &mass_matrix_assembler() const { return mass_matrix_assembler_; }

        inline const std::shared_ptr<Matrix_t> &mass_matrix() const override { return mass_matrix_; }

        bool assemble_mass_matrix() override { return assemble_mass_matrix(*mass_matrix()); }

        void set_mass_matrix(const std::shared_ptr<Matrix_t> &mass_matrix) override { mass_matrix_ = mass_matrix; }

        virtual bool apply(const Vector_t &x, Vector_t &hessian_times_x) const override {
            if (empty(hessian_times_x)) {
                this->space()->create_vector(hessian_times_x);
                rename("hessian_times_x", hessian_times_x);
            } else {
                hessian_times_x *= 0.0;
            }

            bool ok = this->assembler()->apply(x, hessian_times_x);

            if (!ok) return false;

            this->space()->copy_at_constrained_nodes(x, hessian_times_x);

            return ok;
        }

        bool assemble_mass_matrix(Matrix_t &mass_matrix) override {
            this->space()->create_matrix(mass_matrix);

            if (!this->mass_matrix_assembler()->assemble(mass_matrix)) {
                return false;
            }

            rename("mass_matrix", mass_matrix);
            return true;
        }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            ensure_gradient(g);

            assert(this->assembler());

            if (!this->assembler()->assemble(x, g)) {
                return false;
            }

            if (must_apply_constraints_) {
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        bool update(const Vector_t &) override { return true; }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            ensure_hessian(H);

            assert(this->assembler());

            if (!this->assembler()->assemble(x, H)) {
                return false;
            }

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
            }

            return true;
        }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            ensure_hessian(H);
            ensure_gradient(g);

            if (!this->assembler()->assemble(x, H, g)) {
                return false;
            }

            if (export_tensors_) {
                rename("Hu", H);
                write("load_Hu.m", H);

                rename("gu", g);
                write("load_gu.m", g);
            }

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
                this->space()->apply_zero_constraints(g);
            }

            if (export_tensors_) {
                rename("H", H);
                write("load_H.m", H);

                rename("g", g);
                write("load_g.m", g);
            }

            return true;
        }

        inline void create_solution_vector(Vector_t &x) override {
            if (empty(x)) {
                this->space()->create_vector(x);
                x.set(0.0);
            }

            apply_constraints(x);
        }

        inline void apply_constraints(Vector_t &x) const override { this->space()->apply_constraints(x); }

        virtual void set_environment(const std::shared_ptr<Environment_t> &env) override {
            this->assembler()->set_environment(env);
            this->mass_matrix_assembler()->set_environment(env);
        }

        inline const std::shared_ptr<OmniAssembler_t> &assembler() const /*override*/ { return assembler_; }
        inline const std::shared_ptr<FunctionSpace> &space() const override { return space_; }

        bool is_time_dependent() const override { return false; }

        inline bool verbose() const { return verbose_; }
        inline void verbose(const bool val) { verbose_ = val; }

        bool report_solution(const Vector_t &x) override { return space_->write(output_path_, x); }

    protected:
        void ensure_gradient(Vector_t &g) const {
            if (empty(g)) {
                this->space()->create_vector(g);
                rename("g", g);
            } else {
                g *= 0.0;
            }
        }

        void ensure_hessian(Matrix_t &H) const {
            if (empty(H)) {
                this->space()->create_matrix(H);
                rename("H", H);
            } else {
                if (H.is_assembled()) {
                    H *= 0.0;
                }
            }
        }

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<OmniAssembler_t> assembler_;
        std::shared_ptr<OmniAssembler_t> mass_matrix_assembler_;
        std::shared_ptr<Matrix_t> mass_matrix_;

        bool must_apply_constraints_{true};
        bool verbose_{false};
        bool export_tensors_{false};

        utopia::Path output_path_{"output.e"};
    };

    template <class FunctionSpace>
    class TimeDependentFunction : public FEFunctionInterface<FunctionSpace> {
    public:
        using Super = utopia::FEFunctionInterface<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using SimulationTime = utopia::SimulationTime<Scalar_t>;

        using IO_t = utopia::IO<FunctionSpace>;

        TimeDependentFunction(const std::shared_ptr<FunctionSpace> &space)
            : fe_function_(std::make_shared<FEModelFunction_t>(space)) {
            fe_function_->must_apply_constraints_to_assembled(false);
        }

        TimeDependentFunction(const std::shared_ptr<FEFunctionInterface_t> &fe_function) : fe_function_(fe_function) {
            fe_function_->must_apply_constraints_to_assembled(false);
        }

        virtual ~TimeDependentFunction() = default;

        void set_environment(const std::shared_ptr<Environment_t> &env) override { fe_function_->set_environment(env); }

        bool update_IVP(const Vector_t &x) override {
            time()->update();
            space()->update(*time());

            if (fe_function_) fe_function_->update_IVP(x);
            return true;
        }

        void post_solve(Vector_t &x) override { fe_function_->post_solve(x); }

        bool setup_IVP(Vector_t &x) override = 0;

        bool setup_IVP(IO<FunctionSpace> &) override {
            Utopia::Abort("setup_IVP(IO<FunctionSpace> &), not implemented!");
            return false;
        }

        bool is_IVP_solved() override { return time()->finished(); }

        virtual void integrate_gradient(const Vector_t &x, Vector_t &g) const = 0;
        virtual void integrate_hessian(const Vector_t &x, Matrix_t &H) const = 0;

        bool value(const Vector_t &x, Scalar_t &v) const override { return fe_function_->value(x, v); }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            if (!fe_function_->gradient(x, g)) {
                return false;
            }

            if (export_tensors_) {
                rename("g", g);
                write("load_g.m", g);
            }

            integrate_gradient(x, g);

            if (export_tensors_) {
                write("load_gt.m", g);
            }

            if (must_apply_constraints_) {
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        inline void create_solution_vector(Vector_t &x) override { fe_function_->create_solution_vector(x); }
        inline void apply_constraints(Vector_t &x) const override { fe_function_->apply_constraints(x); }

        bool update(const Vector_t &x) override { return fe_function_->update(x); }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            if (!fe_function_->hessian(x, H)) {
                return false;
            }

            if (export_tensors_) {
                rename("H", H);
                write("load_H.m", H);
            }

            integrate_hessian(x, H);

            if (export_tensors_) {
                write("load_Ht.m", H);
            }

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
            }

            return true;
        }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            if (!fe_function_->hessian_and_gradient(x, H, g)) {
                return false;
            }

            if (export_tensors_) {
                rename("H", H);
                rename("g", g);
                write("load_H.m", H);
                write("load_g.m", g);
            }

            integrate_gradient(x, g);
            integrate_hessian(x, H);

            if (export_tensors_) {
                write("load_Ht.m", H);
                write("load_gt.m", g);
            }

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        inline const std::shared_ptr<Matrix_t> &mass_matrix() const override { return fe_function_->mass_matrix(); }

        bool assemble_mass_matrix() override { return fe_function_->assemble_mass_matrix(); }

        void set_mass_matrix(const std::shared_ptr<Matrix_t> &mass_matrix) override {
            fe_function_->set_mass_matrix(mass_matrix);
        }

        bool assemble_mass_matrix(Matrix_t &mass_matrix) override {
            return fe_function_->assemble_mass_matrix(mass_matrix);
        }

        void read(Input &in) override {
            Super::read(in);
            fe_function_->read(in);

            if (!time_) {
                time_ = std::make_shared<SimulationTime>();
            }

            fe_function_->set_time(time_);

            in.require("time", *time_);
            in.get("export_tensors", export_tensors_);

            auto space_name = space()->name();

            if (!space_name.empty()) {
                output_path_ = space_name + ".e";
            }

            in.get("output_path", output_path_);

            in.get("output", [&](Input &node) {
                std::string output_mode = "OVERWRITE";
                node.get("mode", output_mode_);
                node.require("path", output_path_);
            });
        }

        inline Scalar_t delta_time() const { return time()->delta(); }

        inline std::shared_ptr<SimulationTime> &time() {
            assert(time_);
            return time_;
        }

        inline std::shared_ptr<SimulationTime> time() const override {
            assert(time_);
            return time_;
        }

        bool is_time_dependent() const override { return true; }

        // inline const std::shared_ptr<OmniAssembler_t> &assembler() const override { return fe_function_->assembler();
        // }
        inline const std::shared_ptr<FunctionSpace> &space() const override { return fe_function_->space(); }

        inline const std::shared_ptr<FEFunctionInterface_t> function() const { return fe_function_; }
        inline void set_function(const std::shared_ptr<FEFunctionInterface_t> &function) { fe_function_ = function; }

        bool is_linear() const override { return function()->is_linear(); }

        bool is_operator() const override { return function()->is_operator(); }

        Size size() const override { return function()->size(); }

        Size local_size() const override { return function()->local_size(); }

        inline Communicator_t &comm() override { return function()->comm(); }
        inline const Communicator_t &comm() const override { return function()->comm(); }

        inline void must_apply_constraints_to_assembled(const bool val) override { must_apply_constraints_ = val; }

        bool report_solution(const Vector_t &x) override {
            if (!io_) {
                io_ = std::make_shared<IO_t>(*this->space());
                io_->set_output_path(output_path_);
                io_->set_output_mode(output_mode_);
                this->register_output(*io_);
            }

            this->update_output(*io_);

            return io_->write(x, this->time()->step(), this->time()->get());
        }

        virtual const Vector_t &solution() const = 0;

        void set_output_db(const std::shared_ptr<IO_t> &io) override {
            io_ = io;
            this->register_output(*io_);
        }

        void set_time(const std::shared_ptr<SimulationTime> &time) override {
            time_ = time;
            fe_function_->set_time(time_);
        }

    protected:
        inline bool export_tensors() const { return export_tensors_; }
        inline bool must_apply_constraints_to_assembled() const { return must_apply_constraints_; }

    private:
        std::shared_ptr<FEFunctionInterface_t> fe_function_;
        std::shared_ptr<IO_t> io_;

        std::shared_ptr<SimulationTime> time_;

        bool must_apply_constraints_{true};
        bool export_tensors_{false};
        utopia::Path output_path_{"output.e"};
        std::string output_mode_;
    };

}  // namespace utopia

#endif  // UTOPIA_FE_MODEL_FUNCTION_HPP
