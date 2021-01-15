#ifndef UTOPIA_BOX_KORNHUBER_TRUNCATION_HPP
#define UTOPIA_BOX_KORNHUBER_TRUNCATION_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeprecatedHeaders.hpp"
#include "utopia_Function.hpp"
#include "utopia_IPRTruncatedTransfer.hpp"
#include "utopia_IdentityTransfer.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_MultiLevelVariableBoundInterface.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>

namespace utopia {

    template <class Matrix, class Vector>
    class BoxKornhuberTruncation
        : public MultilevelVariableBoundSolverInterface<Matrix, Vector, BoxKornhuberTruncation<Matrix, Vector>> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, BoxKornhuberTruncation<Matrix, Vector>>
            Base;

        BoxKornhuberTruncation(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> &transfer)
            : Base(transfer) {}

        void init_memory_impl(const std::vector<Layout> &layouts) {
            constraints_memory_.init_memory(layouts);
            const SizeType finest_level = layouts.size();

            const SizeType n_levels = layouts.size();
            help_loc_.resize(n_levels);

            for (SizeType l = 0; l < n_levels; l++) {
                help_loc_[l].zeros(layouts[l]);
            }

            if (this->box_constraints_.has_lower_bound()) {
                constraints_memory_.active_lower[finest_level - 1] = *(this->box_constraints_.lower_bound());
            }

            if (this->box_constraints_.has_upper_bound()) {
                constraints_memory_.active_upper[finest_level - 1] = *(this->box_constraints_.upper_bound());
            }
        }

        void init_level_impl(const SizeType &level,
                             const Vector &x_finer_level,
                             const Vector &x_level,
                             const Scalar & /*delta_fine*/) {
            const auto finer_level = level + 1;
            // const auto n_levels = help_loc_.size();

            if (auto *mat_transfer = dynamic_cast<MatrixTransfer<Matrix, Vector> *>(this->transfer_[level].get())) {
                this->help_[finer_level] = constraints_memory_.active_lower[finer_level] - x_finer_level;

                this->help_loc_[finer_level] = constraints_memory_.active_upper[finer_level] - x_finer_level;

                Scalar lb_max = max(this->help_[finer_level]);
                Scalar ub_min = min(this->help_loc_[finer_level]);

                {
                    // TODO(Alena) can this be done with a parallel for?

                    Read<Matrix> r_A(mat_transfer->R());
                    Write<Vector> rw_c_l(constraints_memory_.active_lower[level]);
                    Read<Vector> rw_f_l(this->help_[finer_level]);

                    Write<Vector> rw_c_u(constraints_memory_.active_upper[level]);
                    Read<Vector> rw_f_u(this->help_loc_[finer_level]);

                    Range rr = range(constraints_memory_.active_lower[level]);
                    Range rr_fine_level = range(this->help_[finer_level]);

                    for (auto i = rr.begin(); i != rr.end(); ++i) {
                        RowView<const Matrix> row_view(mat_transfer->R(), i);
                        decltype(i) n_values = row_view.n_values();

                        Scalar max_value = -9e20;
                        Scalar min_value = 9e20;
                        for (auto index = 0; index < n_values; ++index) {
                            const decltype(i) j = row_view.col(index);
                            if (rr_fine_level.inside(j)) {
                                // std::cout<<"----- yes, in...\n";
                                Scalar val_cons_fine_lb = this->help_[finer_level].get(j);

                                // test if node active =>  truncated index, do not take into
                                // account
                                if (device::abs(val_cons_fine_lb) < 1e-14) {
                                    val_cons_fine_lb = -9e9;
                                }

                                max_value = (max_value > val_cons_fine_lb) ? max_value : val_cons_fine_lb;

                                Scalar val_cons_fine_ub = this->help_loc_[finer_level].get(j);

                                // test if node active =>  truncated index, do not take into
                                // account
                                if (device::abs(val_cons_fine_ub) < 1e-14) {
                                    val_cons_fine_ub = 9e9;
                                }

                                min_value = (min_value < val_cons_fine_ub) ? min_value : val_cons_fine_ub;

                            } else {
                                // std::cout<<"----- yes, out...\n";
                                max_value = (max_value > lb_max) ? max_value : lb_max;
                                min_value = (min_value < ub_min) ? min_value : ub_min;
                            }
                        }

                        constraints_memory_.active_lower[level].set(i, max_value);
                        constraints_memory_.active_upper[level].set(i, min_value);
                    }
                }  // R/W lock

                constraints_memory_.active_lower[level] += x_level;
                constraints_memory_.active_upper[level] += x_level;

            }  // dynamic cast test
            else {
                utopia_error(
                    "TRBoundsKornhuberTruncation:: Use IPRTruncatedTransfer! "
                    "\n ");
                std::cout << "--------- error ---------- \n";
            }
        }

        const Vector &active_upper(const SizeType &level) { return constraints_memory_.active_upper[level]; }

        const Vector &active_lower(const SizeType &level) { return constraints_memory_.active_lower[level]; }

    private:
        ConstraintsLevelMemory<Vector> constraints_memory_;
        std::vector<Vector> help_loc_;
    };

}  // namespace utopia
#endif  // UTOPIA_BOX_KORNHUBER_TRUNCATION_HPP
