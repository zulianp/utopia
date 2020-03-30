#ifndef UTOPIA_BOX_GELMAN_MANDEL_HPP
#define UTOPIA_BOX_GELMAN_MANDEL_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_IdentityTransfer.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Matrix, class Vector>
    class BoxGelmanMandel : public MultilevelVariableBoundSolverInterface<Matrix, Vector, BoxGelmanMandel<Matrix, Vector> >
    {
        public:
            using Scalar   = typename Traits<Vector>::Scalar;
            using SizeType = typename Traits<Vector>::SizeType;
            using Layout   = typename Traits<Vector>::Layout;

            typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, BoxGelmanMandel<Matrix, Vector> > Base;

            BoxGelmanMandel(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> &transfer) : Base(transfer)
            {}

            void init_memory_impl(const std::vector<Layout> &layouts)
            {
                constraints_memory_.init_memory(layouts);
                const SizeType finest_level = layouts.size();

                if(this->box_constraints_.has_lower_bound()){
                    constraints_memory_.active_lower[finest_level-1] = *(this->box_constraints_.lower_bound());
                }

                if(this->box_constraints_.has_upper_bound()){
                    constraints_memory_.active_upper[finest_level-1] = *(this->box_constraints_.upper_bound());
                }
            }

            void init_level_impl(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine)
            {
                auto finer_level = level + 1;
                Scalar I_inf_norm = this->transfer_[level]->projection_inf_norm();

                //////////////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////// lower bound /////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////////////
                this->help_[finer_level] = constraints_memory_.active_lower[finer_level] - x_finer_level;
                Scalar max_val = max(this->help_[finer_level]);
                Scalar lower_multiplier = 1.0/I_inf_norm * max_val;

                {
                    auto d_x      = const_device_view(x_level);
                    parallel_each_write(constraints_memory_.active_lower[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        return d_x.get(i) + lower_multiplier;
                    });
                }

                //////////////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////// upper bound /////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////////////
                this->help_[finer_level] = constraints_memory_.active_upper[finer_level] - x_finer_level;
                Scalar min_val = min(this->help_[finer_level]);
                Scalar upper_multiplier = 1.0/I_inf_norm * min_val;
                {
                    auto d_x      = const_device_view(x_level);
                    parallel_each_write(constraints_memory_.active_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        return d_x.get(i) + upper_multiplier;
                    });
                }

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
    };


}
#endif //UTOPIA_BOX_GELMAN_MANDEL_HPP
