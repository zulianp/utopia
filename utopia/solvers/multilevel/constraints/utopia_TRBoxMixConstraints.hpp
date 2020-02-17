#ifndef UTOPIA_BOX_TR_CONSTRAINTS_COMBINED_HPP
#define UTOPIA_BOX_TR_CONSTRAINTS_COMBINED_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_Algorithms.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Matrix, class Vector, class TRConstraints, class BoxConstraints>
    class TRBoxMixConstraints : public MultilevelVariableBoundSolverInterface<Matrix, Vector, TRBoxMixConstraints<Matrix, Vector, TRConstraints, BoxConstraints> >
    {
        public:
            typedef UTOPIA_SCALAR(Vector)                           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

            typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, TRBoxMixConstraints<Matrix, Vector, TRConstraints, BoxConstraints> > Base; 

            TRBoxMixConstraints(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> & transfer) : Base(transfer), 
                                                                                                            tr_bounds_(transfer), 
                                                                                                            box_bounds_(transfer)
            {

            }

            void init_memory_impl(const std::vector<SizeType> & n_dofs_)
            {
                constraints_memory_.init_memory(n_dofs_); 
                tr_bounds_.init_memory(n_dofs_); 
                box_bounds_.init_memory(n_dofs_); 

                const SizeType finest_level = n_dofs_.size(); 
                if(this->box_constraints_.has_lower_bound()){
                    constraints_memory_.active_lower[finest_level] = *(this->box_constraints_.lower_bound());
                }

                if(this->box_constraints_.has_upper_bound()){
                    constraints_memory_.active_upper[finest_level] = *(this->box_constraints_.upper_bound());
                }                   
            }

            void init_level_impl(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine)
            {
                tr_bounds_.init_level(level, x_finer_level, x_level, delta_fine);            
                box_bounds_.init_level(level, x_finer_level, x_level, delta_fine); 

                // intersect  lower bounds
                {
                    auto d_tr_lower      = const_device_view(tr_bounds_.active_lower(level));
                    auto d_box_lower     = const_device_view(box_bounds_.active_lower(level));

                    parallel_each_write(constraints_memory_.active_lower[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        return device::max(d_tr_lower.get(i), d_box_lower.get(i)); 
                    });   
                }  

                // intersect  upper bounds
                {
                    auto d_tr_upper      = const_device_view(tr_bounds_.active_upper(level));
                    auto d_box_upper     = const_device_view(box_bounds_.active_upper(level));

                    parallel_each_write(constraints_memory_.active_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        return device::min(d_tr_upper.get(i), d_box_upper.get(i)); 
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
            TRConstraints tr_bounds_; 
            BoxConstraints box_bounds_; 
    };


    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
    template<class Matrix, class Vector>
    using TRGrattonBoxGelmanMandel = utopia::TRBoxMixConstraints<Matrix, Vector, TRBoundsGratton<Matrix, Vector>, BoxGelmanMandel<Matrix, Vector> >;

    template<class Matrix, class Vector>
    using TRGrattonBoxKornhuber = utopia::TRBoxMixConstraints<Matrix, Vector, TRBoundsGratton<Matrix, Vector>, BoxKornhuber<Matrix, Vector> >;

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    template<class Matrix, class Vector>
    using TRKornhuberBoxGelmanMandel = utopia::TRBoxMixConstraints<Matrix, Vector, TRBoundsKornhuber<Matrix, Vector>, BoxGelmanMandel<Matrix, Vector> >;

    template<class Matrix, class Vector>
    using TRKornhuberBoxKornhuber = utopia::TRBoxMixConstraints<Matrix, Vector, TRBoundsKornhuber<Matrix, Vector>, BoxKornhuber<Matrix, Vector> >;

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
    template<class Matrix, class Vector>
    using TRGelmanMandelBoxGelmanMandel = utopia::TRBoxMixConstraints<Matrix, Vector, TRBoundsGelmanMandel<Matrix, Vector>, BoxGelmanMandel<Matrix, Vector> >;

    template<class Matrix, class Vector>
    using TRGelmanMandelBoxKornhuber = utopia::TRBoxMixConstraints<Matrix, Vector, TRBoundsGelmanMandel<Matrix, Vector>, BoxKornhuber<Matrix, Vector> >;    



}
#endif //UTOPIA_BOX_TR_CONSTRAINTS_COMBINED_HPP
