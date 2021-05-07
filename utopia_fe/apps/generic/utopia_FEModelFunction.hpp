#ifndef UTOPIA_FE_MODEL_FUNCTION_HPP
#define UTOPIA_FE_MODEL_FUNCTION_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_MatrixTransformer.hpp"
#include "utopia_ProblemBase.hpp"

namespace utopia {

    template <class FunctionSpace>
    class FEModelFunction
        : public Function<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;

        FEModelFunction(const std::shared_ptr<FunctionSpace> &space)
            : space_(space), assembler_(std::make_shared<OmniAssembler_t>(space)) {}

        virtual ~FEModelFunction() = default;

        void read(Input &in) override {
            in.get("verbose", verbose_);
            in.get("assembly", *assembler_);
        }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            ensure_gradient(g);

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

        inline void create_solution_vector(Vector_t &x) {
            if (empty(x)) {
                this->space()->create_vector(x);
                x.set(0.0);
            }

            apply_constraints(x);
        }

        inline void apply_constraints(Vector_t &x) const { this->space()->apply_constraints(x); }

        virtual void set_environment(const std::shared_ptr<Environment_t> &env) {
            this->assembler()->set_environment(env);
        }

        inline const std::shared_ptr<OmniAssembler_t> &assembler() const { return assembler_; }
        inline const std::shared_ptr<FunctionSpace> &space() const { return space_; }

        inline bool verbose() const { return verbose_; }
        inline void verbose(const bool val) const { verbose_ = val; }

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

        inline void must_apply_constraints_to_assembled(const bool val) { must_apply_constraints_ = val; }

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<OmniAssembler_t> assembler_;
        bool verbose_{false};
        bool must_apply_constraints_{true};
    };

    template <class FunctionSpace>
    class TimeDependentFunction : public FEModelFunction<FunctionSpace> {
    public:
        using Super = utopia::FEModelFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using Environment_t = utopia::Environment<FunctionSpace>;

        TimeDependentFunction(const std::shared_ptr<FunctionSpace> &space)
            : Super(space),
              mass_matrix_assembler_(std::make_shared<OmniAssembler_t>(space)),
              mass_matrix_(std::make_shared<Matrix_t>()) {
            this->must_apply_constraints_to_assembled(false);
        }

        virtual ~TimeDependentFunction() = default;

        void set_environment(const std::shared_ptr<Environment_t> &env) override {
            Super::set_environment(env);
            this->mass_matrix_assembler()->set_environment(env);
        }

        virtual bool update_IVP(const Vector_t &x) = 0;
        virtual bool setup_IVP(Vector_t &x) = 0;

        void read(Input &in) override {
            Super::read(in);
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

        bool assemble_mass_matrix() {
            this->space()->create_matrix(*mass_matrix());

            if (!this->mass_matrix_assembler()->assemble(*mass_matrix())) {
                return false;
            }

            rename("mass_matrix", *mass_matrix());
            return true;
        }

    private:
        std::shared_ptr<OmniAssembler_t> mass_matrix_assembler_;
        std::shared_ptr<Matrix_t> mass_matrix_;
        Scalar_t delta_time_{0.1};
    };

}  // namespace utopia

#endif  // UTOPIA_FE_MODEL_FUNCTION_HPP
