#ifndef UTOPIA_IDENTITY_CONSTRAINTS_HPP
#define UTOPIA_IDENTITY_CONSTRAINTS_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_Algorithms.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Matrix, class Vector>
    class IdentityConstraints : public MultilevelVariableBoundSolverInterface<Matrix, Vector, IdentityConstraints<Matrix, Vector> >
    {
        public:
            using Scalar   = typename Traits<Vector>::Scalar;
            using SizeType = typename Traits<Vector>::SizeType;
            using Layout   = typename Traits<Vector>::Layout;


            typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, IdentityConstraints<Matrix, Vector> > Base;

            IdentityConstraints(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> & transfer) : Base(transfer)
            {
                // TODO:: move to checck innitialization
                for (auto l=0; l < transfer.size(); l++){
                    if (auto* id_transfer = dynamic_cast<IdentityTransfer<Matrix, Vector>*>(transfer[l].get())) {
                        // utopia_error("IdentityConstraints, termination due to incorrect setup. ");
                    } else {
                        utopia_error("IdentityConstraints, termination due to incorrect setup. ");
                    }
                }
            }


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

            void init_level_impl(const SizeType& level,
                                 const Vector& x_finer_level,
                                 const Vector& /*x_level*/,
                                 const Scalar& delta_fine) {
                auto finer_level = level + 1;

                {
                    auto d_x_finer      = const_device_view(x_finer_level);
                    auto d_tr_lb        = const_device_view(constraints_memory_.active_lower[finer_level]);
                    auto d_tr_ub        = const_device_view(constraints_memory_.active_upper[finer_level]);

                    parallel_each_write(constraints_memory_.active_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        auto val1 = d_x_finer.get(i) + delta_fine;
                        auto val2 = d_tr_ub.get(i);

                        return device::min(val1, val2);
                    });

                    parallel_each_write(constraints_memory_.active_lower[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        auto val1 = d_x_finer.get(i) - delta_fine;
                        auto val2 = d_tr_lb.get(i);

                        return device::max(val1, val2);
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
#endif //UTOPIA_IDENTITY_CONSTRAINTS_HPP
