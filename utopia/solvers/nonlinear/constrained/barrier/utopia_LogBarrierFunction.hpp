#ifndef UTOPIA_LOG_BARRIER_FUNCTION_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_HPP

#include "utopia_make_unique.hpp"

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include "utopia_LogBarrierFunctionBase.hpp"

#include <limits>

namespace utopia {

    template <class Matrix, class Vector>
    class LogBarrier : public LogBarrierBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = utopia::LogBarrierBase<Matrix, Vector>;

        LogBarrier() = default;
        explicit LogBarrier(const std::shared_ptr<BoxConstraints> &box) : Super(box) {}

        bool barrier_hessian(const Vector &x, Vector &H) const override {
            Vector diff;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                diff = pow2(diff);
                H += 1. / diff;
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                diff = pow2(diff);
                H += 1. / diff;
            }

            return true;
        }

        bool barrier_gradient(const Vector &x, Vector &g) const override {
            Vector diff;
            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                g += 1. / diff;
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                g -= 1. / diff;
            }

            return true;
        }

        bool barrier_value(const Vector &x, Vector &value) const override {
            if (this->box_->has_upper_bound()) {
                value -= logn(*this->box_->upper_bound() - x);
            }

            if (this->box_->has_lower_bound()) {
                value += logn(x - *this->box_->lower_bound());
            }

            return true;
        }

        inline std::string function_type() const override { return "LogBarrier"; }

        void read(Input &in) override { Super::read(in); }
    };

    template <class Matrix, class Vector>
    class LogBarrierFunction : public Function<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = Function;
        using LogBarrierBase = utopia::LogBarrierBase<Matrix, Vector>;
        using Transformation = utopia::Transformation<Matrix, Vector>;

        void set_barrier(const std::shared_ptr<LogBarrierBase> &barrier) { this->barrier_ = barrier; }

        const std::shared_ptr<LogBarrierBase> &barrier() const {
            assert(barrier_);
            return this->barrier_;
        }

        void set_transform(const std::shared_ptr<Transformation> &t) { barrier()->set_transform(t); }

        void set_scaling_matrix(const std::shared_ptr<Matrix> &scaling_matrix) {
            barrier()->set_scaling_matrix(scaling_matrix);
        }

        virtual bool extend_hessian_diag(const Vector &x, Vector &h) const {
            if (barrier()) {
                return barrier()->hessian(x, h);
            }
            return false;
        }

        virtual bool extend_hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const {
            if (barrier()) {
                return barrier()->hessian_and_gradient(x, H, g);
            }
            return false;
        }

        virtual bool extend_hessian(const Vector &x, Matrix &H) const {
            if (barrier()) {
                return barrier()->hessian(x, H);
            }
            return false;
        }

        virtual bool extend_gradient(const Vector &x, Vector &g) const {
            if (barrier()) {
                return barrier()->gradient(x, g);
            }
            return false;
        }

        virtual bool extend_value(const Vector &x, Scalar &value) const {
            if (barrier()) {
                return barrier()->value(x, value);
            }
            return false;
        }

        virtual bool extend_project_onto_feasibile_region(Vector &x) const {
            if (barrier()) return barrier()->project_onto_feasibile_region(x);
            return false;
        }

        virtual std::string function_type() const { return barrier()->function_type(); }

        void read(Input &in) override {
            if (barrier()) barrier()->read(in);
            reset();
        }

        inline void set_unconstrained_function(const std::shared_ptr<Function> &unconstrained) {
            unconstrained_ = unconstrained;
        }

        void set_box_constraints(const std::shared_ptr<BoxConstraints> &box) { barrier()->set_box_constraints(box); }
        void set_selection(const std::shared_ptr<Vector> &boolean_selector) {
            barrier()->set_selection(boolean_selector);
        }

        bool hessian_diag(const Vector &x, Vector &h) const {
            if (unconstrained_) {
                assert(false);
                Utopia::Abort("IMPLEMENT ME");
            }

            return barrier()->hessian(x, h);
        }

        bool hessian(const Vector &x, Matrix &H) const final {
            if (unconstrained_) {
                if (!unconstrained_->hessian(x, H)) {
                    return false;
                }
            }

            return extend_hessian(x, H);
        }

        bool hessian(const Vector &x, Matrix &H, Matrix &preconditioner) const final {
            if (unconstrained_) {
                if (!unconstrained_->hessian(x, H, preconditioner)) {
                    return false;
                }
            }

            return extend_hessian(x, H);
        }

        bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const final {
            if (unconstrained_) {
                if (!unconstrained_->hessian_and_gradient(x, H, g)) {
                    return false;
                }
            }

            return extend_hessian_and_gradient(x, H, g);
        }

        bool project_onto_feasibile_region(Vector &x) const final {
            return barrier()->project_onto_feasibile_region(x);
        }

        bool hessian_and_gradient(const Vector &x, Matrix &H, Matrix &preconditioner, Vector &g) const override {
            if (unconstrained_) {
                if (!unconstrained_->hessian_and_gradient(x, H, preconditioner, g)) {
                    return false;
                }
            }

            return extend_hessian_and_gradient(x, H, g);
        }

        inline bool has_preconditioner() const override { return unconstrained_->has_preconditioner(); }

        bool initialize_hessian(Matrix &H, Matrix &preconditioner) const override {
            if (unconstrained_) {
                if (!unconstrained_->initialize_hessian(H, preconditioner)) {
                    return false;
                }
            }

            return true;
        }

        bool value(const Vector &x, Scalar &value) const override {
            if (unconstrained_) {
                if (!unconstrained_->value(x, value)) {
                    value = 0.0;
                }
            }

            return extend_value(x, value);
        }

        bool gradient(const Vector &x, Vector &g) const override {
            if (unconstrained_) {
                if (!unconstrained_->gradient(x, g)) {
                    return false;
                }
            }

            return extend_gradient(x, g);
        }

        bool update(const Vector &x) override {
            if (unconstrained_) {
                if (!unconstrained_->update(x)) {
                    return false;
                }
            }

            update_barrier();
            return true;
        }

        inline void set_barrier_parameter(const Scalar value) {
            barrier()->barrier_parameter_ = value;
            barrier()->current_barrier_parameter_ = value;
        }

        inline void set_barrier_parameter_shrinking_factor(const Scalar value) {
            barrier()->barrier_parameter_shrinking_factor_ = value;
        }

        inline void set_min_barrier_parameter(const Scalar value) { barrier()->min_barrier_parameter_ = value; }

        virtual void reset() {
            if (!barrier_) {
                barrier_ = std::make_shared<LogBarrier<Matrix, Vector>>();
            }

            if (barrier()) barrier()->reset();
        }

        inline bool verbose() const {
            if (barrier())
                return barrier()->verbose();
            else
                return false;
        }

        inline bool has_transformation() const { return barrier()->has_transformation(); }

        void update_barrier() {
            if (barrier()) barrier()->update_barrier();
        }

        void compute_diff_upper_bound(const Vector &x, Vector &diff) const {
            if (barrier()) barrier()->compute_diff_upper_bound(x, diff);
        }

        void compute_diff_lower_bound(const Vector &x, Vector &diff) const {
            if (barrier()) barrier()->compute_diff_lower_bound(x, diff);
        }

        void auto_selector(const bool val) { barrier()->auto_selector(val); }

        LogBarrierFunction() : barrier_(std::make_shared<LogBarrier<Matrix, Vector>>()) {}

        LogBarrierFunction(const std::shared_ptr<Function> &unconstrained,
                           const std::shared_ptr<LogBarrierBase> &barrier)
            : unconstrained_(unconstrained), barrier_(barrier) {}

        LogBarrierFunction(const std::shared_ptr<LogBarrierBase> &barrier) : barrier_(barrier) {}

        UTOPIA_NVCC_PRIVATE
        std::shared_ptr<Function> unconstrained_;
        std::shared_ptr<LogBarrierBase> barrier_;
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_HPP
