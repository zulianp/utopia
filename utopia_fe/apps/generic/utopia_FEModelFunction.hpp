#ifndef UTOPIA_FE_MODEL_FUNCTION_HPP
#define UTOPIA_FE_MODEL_FUNCTION_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_MatrixTransformer.hpp"
#include "utopia_ProblemBase.hpp"

#include <limits>

namespace utopia {

    template <class FunctionSpace>
    class FEFunctionInterface
        : public Function<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;

        virtual ~FEFunctionInterface() = default;

        virtual void create_solution_vector(Vector_t &x) = 0;
        virtual void apply_constraints(Vector_t &x) const = 0;
        virtual void set_environment(const std::shared_ptr<Environment_t> &env) = 0;

        virtual const std::shared_ptr<OmniAssembler_t> &assembler() const = 0;
        virtual const std::shared_ptr<FunctionSpace> &space() const = 0;

        bool initialize_hessian(Matrix_t &H, Matrix_t &) const override {
            space()->create_matrix(H);
            return true;
        }

        virtual bool is_time_dependent() const = 0;
        virtual bool is_linear() const = 0;
    };

    template <class FunctionSpace>
    class FEModelFunction final : public FEFunctionInterface<FunctionSpace> {
    public:
        using Super = utopia::FEFunctionInterface<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;

        FEModelFunction(const std::shared_ptr<FunctionSpace> &space)
            : space_(space), assembler_(std::make_shared<OmniAssembler_t>(space)) {}

        virtual ~FEModelFunction() = default;

        bool is_linear() const override { return assembler()->is_linear(); }

        void read(Input &in) override {
            Super::read(in);
            in.get("assembly", *assembler_);
            in.get("verbose", verbose_);
        }

        bool value(const Vector_t & /*point*/, Scalar_t &value) const override {
            // assert(false && "IMPLEMENT ME");
            value = std::numeric_limits<Scalar_t>::signaling_NaN();
            return false;
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

        bool update(const Vector_t &x) override { return true; }

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

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
                this->space()->apply_zero_constraints(g);
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
        }

        inline const std::shared_ptr<OmniAssembler_t> &assembler() const override { return assembler_; }
        inline const std::shared_ptr<FunctionSpace> &space() const override { return space_; }

        bool is_time_dependent() const override { return false; }

        inline void must_apply_constraints_to_assembled(const bool val) { must_apply_constraints_ = val; }

        inline bool verbose() const { return verbose_; }
        inline void verbose(const bool val) const { verbose_ = val; }

    private:
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

        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<OmniAssembler_t> assembler_;
        bool must_apply_constraints_{true};
        bool verbose_{false};
    };

    template <class FunctionSpace>
    class TimeDependentFunction : public FEFunctionInterface<FunctionSpace> {
    public:
        using Super = utopia::FEFunctionInterface<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        ;

        TimeDependentFunction(const std::shared_ptr<FunctionSpace> &space)
            : fe_function_(std::make_shared<FEModelFunction_t>(space)),
              mass_matrix_assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_(std::make_shared<Matrix_t>()) {
            fe_function_->must_apply_constraints_to_assembled(false);
        }

        virtual ~TimeDependentFunction() = default;

        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            fe_function_->set_environment(env);
            this->mass_matrix_assembler()->set_environment(env);
        }

        virtual bool update_IVP(const Vector_t &x) = 0;
        virtual bool setup_IVP(Vector_t &x) = 0;

        virtual void integrate_gradient(const Vector_t &x, Vector_t &g) const = 0;
        virtual void integrate_hessian(const Vector_t &x, Matrix_t &H) const = 0;

        bool value(const Vector_t &x, Scalar_t &v) const override { return fe_function_->value(x, v); }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            if (!this->assembler()->assemble(x, g)) {
                return false;
            }

            integrate_gradient(x, g);

            if (must_apply_constraints_) {
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        inline void create_solution_vector(Vector_t &x) override { fe_function_->create_solution_vector(x); }
        inline void apply_constraints(Vector_t &x) const override { fe_function_->apply_constraints(x); }

        bool update(const Vector_t &x) override { return fe_function_->update(x); }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            if (!this->assembler()->assemble(x, H)) {
                return false;
            }

            integrate_hessian(x, H);

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
            }

            return true;
        }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            if (!this->assembler()->assemble(x, H, g)) {
                return false;
            }

            integrate_gradient(x, g);
            integrate_hessian(x, H);

            if (must_apply_constraints_) {
                this->space()->apply_constraints(H);
                this->space()->apply_zero_constraints(g);
            }

            return true;
        }

        void read(Input &in) override {
            Super::read(in);
            fe_function_->read(in);

            in.get("delta_time", delta_time_);

            bool user_defined_mass = false;
            in.get("mass", [this, &user_defined_mass](Input &node) {
                this->mass_matrix_assembler()->read(node);
                user_defined_mass = true;
            });

            if (!user_defined_mass) {
                auto params = param_list(param("material",
                                               param_list(param("type", "Mass"),                         //
                                                          param("lumped", true),                         //
                                                          param("n_components", this->space()->n_var())  //
                                                          )));
                this->mass_matrix_assembler()->read(params);
            }
        }

        inline Scalar_t delta_time() const { return delta_time_; }
        inline const std::shared_ptr<Matrix_t> &mass_matrix() const { return mass_matrix_; }
        inline const std::shared_ptr<OmniAssembler_t> &mass_matrix_assembler() const { return mass_matrix_assembler_; }

        bool is_time_dependent() const override { return true; }

        bool assemble_mass_matrix() {
            this->space()->create_matrix(*mass_matrix());

            if (!this->mass_matrix_assembler()->assemble(*mass_matrix())) {
                return false;
            }

            rename("mass_matrix", *mass_matrix());
            return true;
        }

        inline const std::shared_ptr<OmniAssembler_t> &assembler() const override { return fe_function_->assembler(); }
        inline const std::shared_ptr<FunctionSpace> &space() const override { return fe_function_->space(); }

        inline const std::shared_ptr<FEModelFunction_t> function() const { return fe_function_; }

        bool is_linear() const override { return function()->is_linear(); }

    protected:
        inline void must_apply_constraints_to_assembled(const bool val) { must_apply_constraints_ = val; }

    private:
        std::shared_ptr<FEModelFunction_t> fe_function_;
        std::shared_ptr<OmniAssembler_t> mass_matrix_assembler_;
        std::shared_ptr<Matrix_t> mass_matrix_;
        Scalar_t delta_time_{0.1};
        bool must_apply_constraints_{true};
    };

}  // namespace utopia

#endif  // UTOPIA_FE_MODEL_FUNCTION_HPP