#ifndef UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
#define UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_Function.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>

namespace utopia {
    template <class Vector>
    class VariableBoundSolverInterface : virtual public MemoryInterface<Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using BoxConstraints = utopia::BoxConstraints<Vector>;

    public:
        VariableBoundSolverInterface() = default;
        ~VariableBoundSolverInterface() override = default;

        VariableBoundSolverInterface(const VariableBoundSolverInterface &other) : constraints_(other.constraints_) {}

        virtual bool set_box_constraints(const BoxConstraints &box) {
            constraints_ = box;
            return true;
        }

        const BoxConstraints &get_box_constraints() const { return constraints_; }
        BoxConstraints &get_box_constraints() { return constraints_; }

        virtual const Vector &get_upper_bound() const {
            if (!constraints_.has_upper_bound()) {
                utopia_error("VariableBoundSolverInterface::upper bound does not exist. \n");
            }

            return *constraints_.upper_bound();
        }

        virtual const Vector &get_lower_bound() const {
            if (!constraints_.has_lower_bound()) {
                utopia_error("VariableBoundSolverInterface::lower bound does not exist. \n");
            }

            return *constraints_.lower_bound();
        }

        virtual std::shared_ptr<Vector> &upper_bound() { return constraints_.upper_bound(); }

        virtual std::shared_ptr<Vector> &lower_bound() { return constraints_.lower_bound(); }

        virtual bool has_bound() const { return constraints_.has_bound(); }

        virtual bool has_lower_bound() const { return constraints_.has_lower_bound(); }

        virtual bool has_upper_bound() const { return constraints_.has_upper_bound(); }

        virtual bool has_empty_bounds() const { return constraints_.has_empty_bounds(); }

        virtual void fill_empty_bounds(const Layout &layout) { return constraints_.fill_empty_bounds(layout); }

    protected:
        virtual Scalar criticality_measure_infty(const Vector &x, const Vector &g) {
            help_ = x - g;

            if (!constraints_.has_upper_bound() || !constraints_.has_lower_bound()) {
                this->fill_empty_bounds(layout(x));
            }

            const auto &ub = *constraints_.upper_bound();
            const auto &lb = *constraints_.lower_bound();

            get_projection(lb, ub, help_);
            help_ -= x;

            return norm2(help_);
        }

        virtual void project_gradient(const Vector &x, Vector &g) {
            g = x - g;

            if (!constraints_.has_upper_bound() || !constraints_.has_lower_bound()) {
                this->fill_empty_bounds(layout(x));
            }

            const auto &ub = *constraints_.upper_bound();
            const auto &lb = *constraints_.lower_bound();

            get_projection(lb, ub, g);
            g -= x;
        }

    public:
        // expose it for CUDA
        UTOPIA_DEPRECATED_MSG(
            "Use void project(const Vector &lb, const Vector &ub, Vector &x) instead, which does not need a secondary "
            "vector.")
        bool get_projection(const Vector &x, const Vector &lb, const Vector &ub, Vector &Pc) const {
            if (empty(Pc) || size(Pc) != size(x)) {
                Pc = 0.0 * x;
            }

            assert(x.size() == lb.size());
            assert(x.size() == ub.size());
            assert(x.comm().size() == Pc.comm().size());

            {
                auto d_lb = const_view_device(lb);
                auto d_ub = const_view_device(ub);
                auto d_x = const_view_device(x);

                parallel_each_write(Pc, UTOPIA_LAMBDA(const SizeType i)->Scalar {
                    const Scalar li = d_lb.get(i);
                    const Scalar ui = d_ub.get(i);
                    const Scalar xi = d_x.get(i);
                    return (li >= xi) ? li : ((ui <= xi) ? ui : xi);
                });
            }

            return true;
        }

        void project(const Vector &lb, const Vector &ub, Vector &x) const {
            assert(!empty(lb));
            assert(!empty(ub));

            assert(x.size() == lb.size());
            assert(x.size() == ub.size());
            assert(x.comm().size() == lb.comm().size());

            {
                auto d_lb = const_local_view_device(lb);
                auto d_ub = const_local_view_device(ub);
                auto d_x = local_view_device(x);

                parallel_for(local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    const Scalar li = d_lb.get(i);
                    const Scalar ui = d_ub.get(i);
                    const Scalar xi = d_x.get(i);
                    d_x.set(i, (li >= xi) ? li : ((ui <= xi) ? ui : xi));
                });
            }
        }

        UTOPIA_DEPRECATED_MSG(
            "Use void project(const Vector &lb, const Vector &ub, Vector &x) instead, which uses local indexing "
            "(faster)")
        bool get_projection(const Vector &lb, const Vector &ub, Vector &x) const {
            {
                auto d_lb = const_view_device(lb);
                auto d_ub = const_view_device(ub);

                parallel_transform(x, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi)->Scalar {
                    const Scalar li = d_lb.get(i);
                    const Scalar ui = d_ub.get(i);
                    return (li >= xi) ? li : ((ui <= xi) ? ui : xi);
                });
            }

            return true;
        }

        void make_iterate_feasible(Vector &x) const {
            if (!constraints_.has_upper_bound() && !constraints_.has_lower_bound()) return;

            if (constraints_.has_upper_bound() && constraints_.has_lower_bound()) {
                const auto &ub = *constraints_.upper_bound();
                const auto &lb = *constraints_.lower_bound();

                {
                    auto d_lb = const_view_device(lb);
                    auto d_ub = const_view_device(ub);

                    parallel_transform(x, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi)->Scalar {
                        const Scalar li = d_lb.get(i);
                        const Scalar ui = d_ub.get(i);
                        return (li >= xi) ? li : ((ui <= xi) ? ui : xi);
                    });
                }

            } else if (constraints_.has_upper_bound() && !constraints_.has_lower_bound()) {
                const auto &ub = *constraints_.upper_bound();

                {
                    auto d_ub = const_view_device(ub);

                    parallel_transform(x, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi)->Scalar {
                        const Scalar ui = d_ub.get(i);
                        return (ui <= xi) ? ui : xi;
                    });
                }
            } else {
                const auto &lb = *constraints_.lower_bound();

                {
                    auto d_lb = const_view_device(lb);

                    parallel_transform(x, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi)->Scalar {
                        const Scalar li = d_lb.get(i);
                        return (li >= xi) ? li : xi;
                    });
                }
            }
        }

        virtual const BoxConstraints &merge_pointwise_constraints_with_uniform_bounds(const Vector &x_k,
                                                                                      const Scalar &lb_uniform,
                                                                                      const Scalar &ub_uniform) {
            correction_constraints_.fill_empty_bounds(layout(x_k));

            if (constraints_.has_upper_bound()) {
                auto &ub_merged = *correction_constraints_.upper_bound();
                ub_merged = *constraints_.upper_bound() - x_k;

                {
                    parallel_transform(ub_merged, UTOPIA_LAMBDA(const SizeType &, const Scalar &ub)->Scalar {
                        return device::min(ub, ub_uniform);
                    });
                }
            } else {
                correction_constraints_.upper_bound()->set(ub_uniform);
            }

            if (constraints_.has_lower_bound()) {
                auto &lb_merged = *correction_constraints_.lower_bound();
                lb_merged = *constraints_.lower_bound() - x_k;

                {
                    parallel_transform(lb_merged, UTOPIA_LAMBDA(const SizeType &, const Scalar &lb)->Scalar {
                        return device::max(lb, lb_uniform);
                    });
                }

            } else {
                correction_constraints_.lower_bound()->set(lb_uniform);
            }

            return correction_constraints_;
        }

        virtual const BoxConstraints &build_correction_constraints(const Vector &x_k) {
            correction_constraints_.fill_empty_bounds(layout(x_k));

            if (constraints_.has_upper_bound()) {
                *correction_constraints_.upper_bound() = *constraints_.upper_bound() - x_k;
            } else {
                correction_constraints_.upper_bound()->set(9e12);
            }

            if (constraints_.has_lower_bound()) {
                *correction_constraints_.lower_bound() = *(constraints_.lower_bound()) - x_k;
            } else {
                correction_constraints_.lower_bound()->set(-9e12);
            }

            return correction_constraints_;
        }

        void init_memory(const Layout &layout) override {
            help_.zeros(layout);
            constraints_.fill_empty_bounds(layout);
            correction_constraints_.fill_empty_bounds(layout);
        }

    protected:
        BoxConstraints constraints_;             // variable bound constraints
        BoxConstraints correction_constraints_;  // constraints needed for correction

        // Vector Pc_;
        // Vector xg_;
        Vector help_;
    };

}  // namespace utopia
#endif  // UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
