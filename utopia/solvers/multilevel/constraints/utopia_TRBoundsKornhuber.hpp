#ifndef UTOPIA_TR_BOUNDS_KORNHUBER_HPP
#define UTOPIA_TR_BOUNDS_KORNHUBER_HPP

#include "utopia_Core.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_IdentityTransfer.hpp"

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

                    // to be deleted 
                    // constraints_memory_.active_lower[finer_level].set(-delta_fine);
                    // constraints_memory_.active_upper[finer_level].set(delta_fine);


                    if(MatrixTransfer<Matrix, Vector>* mat_transfer =  dynamic_cast<MatrixTransfer<Matrix, Vector>* > (this->transfer_[level].get()))
                    {
                        // disp(id_transfer->I(), "I");
                        // std::cout<<"loc_size(I).get(0): "<< local_size(id_transfer->I()).get(0) << "  :  "<< local_size(id_transfer->I()).get(1) << "  \n"; 
                        // disp(mat_transfer->R(), "R");
                        // exit(0);


                        {
                            auto d_x_finer      = const_device_view(x_finer_level);
                            auto d_tr_lb        = const_device_view(constraints_memory_.active_lower[finer_level]);
                            auto d_tr_ub        = const_device_view(constraints_memory_.active_upper[finer_level]);

                            parallel_each_write(this->help_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                            {
                                Scalar xi = d_x_finer.get(i); 
                                auto val = xi - delta_fine;
                                auto lbi = d_tr_lb.get(i); 

                                return std::max(lbi, val) - xi;
                            });   

                            parallel_each_write(this->help_loc_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                            {
                                Scalar xi = d_x_finer.get(i); 
                                auto val = xi + delta_fine;
                                auto ubi = d_tr_ub.get(i); 

                                return std::min(ubi, val) - xi;
                            });                               
                        }

                        Scalar lb_max = max(this->help_[finer_level]); 
                        Scalar ub_min = min(this->help_loc_[finer_level]); 

                        {
                            Read<Matrix> r_A(mat_transfer->R());
                            Write<Vector> rw_c_l(constraints_memory_.active_lower[level]);
                            Read<Vector> rw_f_l(constraints_memory_.active_lower[finer_level]);

                            Write<Vector> rw_c_u(constraints_memory_.active_upper[level]);
                            Read<Vector> rw_f_u(constraints_memory_.active_upper[finer_level]);                            


                            Range rr = range(constraints_memory_.active_lower[level]);  

                            for(auto i = rr.begin(); i != rr.end(); ++i) {
                                RowView<const Matrix> row_view(mat_transfer->R(), i);
                                decltype(i) n_values = row_view.n_values();                            

                                // std::cout<<"n_values: "<< n_values <<"  \n"; 

                                Scalar max_value = -9e9; 
                                Scalar min_value = 9e9; 
                                for(auto index = 0; index < n_values; ++index) {
                                    const decltype(i) j = row_view.col(index);

                                    // std::cout<<"i: "<< i << "  j: "<< j << "  \n"; 
                                    // check if inside, else put lb_max 
                                    Scalar val_cons_fine_lb = constraints_memory_.active_lower[finer_level].get(j); 
                                    max_value = (max_value > val_cons_fine_lb) ? max_value : val_cons_fine_lb; 

                                    Scalar val_cons_fine_ub = constraints_memory_.active_upper[finer_level].get(j); 
                                    min_value = (min_value < val_cons_fine_ub) ? min_value : val_cons_fine_ub;                                     
                                }
                                constraints_memory_.active_lower[level].set(i, max_value); 
                                constraints_memory_.active_upper[level].set(i, min_value); 
                            }

                            // disp(constraints_memory_.active_lower[level], "constraints_memory_.active_lower[level]"); 
                            // disp(constraints_memory_.active_upper[level], "constraints_memory_.active_upper[level]"); 
                            // exit(0);

                        }



                    }
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
            ActiveConstraintsLevelMemory<Vector> constraints_memory_; 
            std::vector<Vector> help_loc_; 
    };

}
#endif //UTOPIA_TR_BOUNDS_KORNHUBER_HPP
