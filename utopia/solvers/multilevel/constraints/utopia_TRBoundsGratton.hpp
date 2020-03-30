#ifndef UTOPIA_TR_BOUNDS_GRATTON_HPP
#define UTOPIA_TR_BOUNDS_GRATTON_HPP

#include "utopia_Core.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_IdentityTransfer.hpp"
#include "utopia_Algorithms.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    template<class Matrix, class Vector>
    class TRBoundsGratton : public MultilevelVariableBoundSolverInterface<Matrix, Vector, TRBoundsGratton<Matrix, Vector> >
    {
        public:
            using Scalar   = typename Traits<Vector>::Scalar;
            using SizeType = typename Traits<Vector>::SizeType;
            using Layout   = typename Traits<Vector>::Layout;

            typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, TRBoundsGratton<Matrix, Vector> > Base;

            TRBoundsGratton(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> & transfer) : Base(transfer)
            {

            }


            void init_memory_impl(const std::vector<Layout> &layouts)
            {
                constraints_memory_.init_memory(layouts);

                const auto n_levels = layouts.size();
                help_loc_.resize(n_levels);
                for(auto l=0; l < n_levels; l++){
                    help_loc_[l].zeros(layouts[l]);
                }

            }

            void init_level_impl(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine)
            {
                auto finer_level = level + 1;
                {
                    auto d_x_finer      = const_device_view(x_finer_level);
                    auto d_tr_lb        = const_device_view(constraints_memory_.active_lower[finer_level]);
                    auto d_tr_ub        = const_device_view(constraints_memory_.active_upper[finer_level]);

                    parallel_each_write(this->help_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        auto val = d_x_finer.get(i) - delta_fine;
                        auto lbi = d_tr_lb.get(i);
                        return device::max(lbi, val);
                    });

                    parallel_each_write(help_loc_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        auto val = d_x_finer.get(i) + delta_fine;
                        auto ubi = d_tr_ub.get(i);
                        return device::min(ubi, val);
                    });
                }

                //------------------------ we should take into account  positive and negative elements projection separatelly -----------------
                this->transfer_[level]->project_down_positive_negative(this->help_[finer_level], help_loc_[finer_level], constraints_memory_.active_lower[level]);
                this->transfer_[level]->project_down_positive_negative(help_loc_[finer_level], this->help_[finer_level], constraints_memory_.active_upper[level]);
            }

            const Vector & active_upper(const SizeType & level)
            {
                return constraints_memory_.active_upper[level];
            }

            const Vector & active_lower(const SizeType & level)
            {
                return constraints_memory_.active_lower[level];
            }

        private:
            ConstraintsLevelMemory<Vector> constraints_memory_;
            std::vector<Vector> help_loc_;
    };

}
#endif //UTOPIA_TR_BOUNDS_GRATTON_HPP
