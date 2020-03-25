#ifndef UTOPIA_TR_BOUNDS_KORNHUBER_HPP
#define UTOPIA_TR_BOUNDS_KORNHUBER_HPP

#include "utopia_Core.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_IdentityTransfer.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_IPTransferNested.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    template<class Matrix, class Vector>
    class TRBoundsKornhuber : public MultilevelVariableBoundSolverInterface<Matrix, Vector, TRBoundsKornhuber<Matrix, Vector> >
    {
        public:
            typedef UTOPIA_SCALAR(Vector)                           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

            typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, TRBoundsKornhuber<Matrix, Vector> > Base; 

            TRBoundsKornhuber(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> & transfer) : Base(transfer)
            {

            }


            void init_memory_impl(const std::vector<SizeType> & n_dofs_)
            {
                constraints_memory_.init_memory(n_dofs_); 
                
                const auto n_levels = n_dofs_.size(); 
                help_loc_.resize(n_levels); 
                for(auto l=0; l < n_levels; l++){
                    help_loc_[l] = local_zeros(n_dofs_[l]); 
                }                                
            }

            void init_level_impl(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine)
            {
                auto finer_level = level + 1; 
                if(finer_level == constraints_memory_.active_lower.size()-1){
                    constraints_memory_.active_lower[level].set(-delta_fine);  
                    constraints_memory_.active_lower[level] += x_level; 

                    constraints_memory_.active_upper[level].set(delta_fine); 
                    constraints_memory_.active_upper[level] += x_level; 
                }
                else
                {
                    // TODO:: check for other transfers 
                    // todo:: this version probably does not work on GPU 
                    if(MatrixTransfer<Matrix, Vector>* mat_transfer =  dynamic_cast<MatrixTransfer<Matrix, Vector>* > (this->transfer_[level].get()))
                    // if(IPTransferNested<Matrix, Vector>* mat_transfer =  dynamic_cast<IPTransferNested<Matrix, Vector>* > (this->transfer_[level].get()))
                    {
                        {
                            auto d_x_finer      = const_device_view(x_finer_level);
                            auto d_tr_lb        = const_device_view(constraints_memory_.active_lower[finer_level]);
                            auto d_tr_ub        = const_device_view(constraints_memory_.active_upper[finer_level]);

                            parallel_each_write(this->help_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                            {
                                Scalar xi = d_x_finer.get(i); 
                                auto val = xi - delta_fine;
                                auto lbi = d_tr_lb.get(i); 

                                return device::max(lbi, val) - xi;
                            });   

                            parallel_each_write(this->help_loc_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                            {
                                Scalar xi = d_x_finer.get(i); 
                                auto val = xi + delta_fine;
                                auto ubi = d_tr_ub.get(i); 

                                return device::min(ubi, val) - xi;
                            });                               
                        }

                        Scalar lb_max = max(this->help_[finer_level]); 
                        Scalar ub_min = min(this->help_loc_[finer_level]); 

                        {
                            Read<Matrix> r_A(mat_transfer->R());
                            Write<Vector> rw_c_l(constraints_memory_.active_lower[level]);
                            Read<Vector> rw_f_l(this->help_[finer_level]);

                            Write<Vector> rw_c_u(constraints_memory_.active_upper[level]);
                            Read<Vector> rw_f_u(this->help_loc_[finer_level]);                            

                            Range rr = range(constraints_memory_.active_lower[level]);  
                            Range rr_fine_level = range(this->help_[finer_level]);  

                            for(auto i = rr.begin(); i != rr.end(); ++i) {
                                RowView<const Matrix> row_view(mat_transfer->R(), i);
                                decltype(i) n_values = row_view.n_values();                            

                                Scalar max_value = -9e20; 
                                Scalar min_value = 9e20; 
                                for(auto index = 0; index < n_values; ++index) {
                                    const decltype(i) j = row_view.col(index);
                                    if(rr_fine_level.inside(j)) {
                                        Scalar val_cons_fine_lb = this->help_[finer_level].get(j); 
                                        max_value = (max_value > val_cons_fine_lb) ? max_value : val_cons_fine_lb; 

                                        Scalar val_cons_fine_ub = this->help_loc_[finer_level].get(j); 
                                        min_value = (min_value < val_cons_fine_ub) ? min_value : val_cons_fine_ub;                                     
                                    }
                                    else{
                                        max_value = (max_value > lb_max) ? max_value : lb_max; 
                                        min_value = (min_value < ub_min) ? min_value : ub_min;                                              
                                    }
                                }
                                constraints_memory_.active_lower[level].set(i, max_value); 
                                constraints_memory_.active_upper[level].set(i, min_value); 
                            }
                        } // R/W lock 

                        constraints_memory_.active_lower[level] += x_level; 
                        constraints_memory_.active_upper[level] += x_level;                         
                        
                    } // dynamic cast test 
                    else{
                        utopia_error("TRBoundsKornhuber:: transfer operators not supported. \n "); 
                        std::cout<<"--------- error ---------- \n"; 
                    }

                } // level check 

            } // fun end 

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
#endif //UTOPIA_TR_BOUNDS_KORNHUBER_HPP
