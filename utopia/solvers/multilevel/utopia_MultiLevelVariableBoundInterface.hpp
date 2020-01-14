#ifndef UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP
#define UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"

#include "utopia_IdentityTransfer.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    template<class Matrix, class Vector, class Child>
    class MultilevelVariableBoundSolverInterface 
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
        typedef utopia::BoxConstraints<Vector>                  BoxConstraints;

        using Device   = typename Traits<Vector>::Device;


    public:
        MultilevelVariableBoundSolverInterface(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> & transfer): 
                                                has_box_constraints_(false), 
                                                transfer_(transfer)
        {

        }

        virtual ~MultilevelVariableBoundSolverInterface() = default;

        virtual void set_box_constraints(BoxConstraints & box)
        {
            box_constraints_ = box;
            has_box_constraints_ = true;
        }

        virtual BoxConstraints & get_box_constraints()
        {
            return box_constraints_;
        }

        const Vector & active_upper(const SizeType & level)
        {
            return static_cast<Child*>(this)->active_upper(level); 
        }

        const Vector & active_lower(const SizeType & level)
        {
            return static_cast<Child*>(this)->active_lower(level); 
        }        


        // TODO:: add checks for proper template if has bounds            


      protected:
        void init_memory(const std::vector<SizeType> & n_dofs_)
        {
            // constraints_memory_.init_memory(n_dofs_);

            // SizeType fine_local_size = n_dofs_.back(); 
            // SizeType n_levels =  n_dofs_.size(); 

            // const SizeType fine_level =n_levels -1;
            // const Scalar inf = std::numeric_limits<Scalar>::infinity();

            // if(has_box_constraints_)
            // {
            //     if(box_constraints_.has_upper_bound()){
            //         constraints_memory_.x_upper[fine_level] = *box_constraints_.upper_bound();
            //     }
            //     else{
            //         constraints_memory_.x_upper[fine_level].set(inf); 
            //     }

            //     if(box_constraints_.has_lower_bound()){
            //         constraints_memory_.x_lower[fine_level] = *box_constraints_.lower_bound();
            //     }
            //     else{
            //         constraints_memory_.x_lower[fine_level].set(-1.0 * inf);
            //     }

            //     constraints_memory_.active_upper[fine_level] = constraints_memory_.x_upper[fine_level];
            //     constraints_memory_.active_lower[fine_level] = constraints_memory_.x_lower[fine_level];
            // }
            // else
            // {
            //     constraints_memory_.active_upper[fine_level].set(inf);
            //     constraints_memory_.active_lower[fine_level].set(-1.0 * inf);                
            // }

            // // inherited tr bound constraints...
            // constraints_memory_.tr_upper[fine_level].set(inf);
            // constraints_memory_.tr_lower[fine_level].set(-1.0 * inf);            

            // TODO:: add lower and upper bound on finest level 
            // precompute norms of prolongation operators needed for projections of constraints...
    
            const auto n_levels = n_dofs_.size();             
            help_.resize(n_levels); 

            for(auto l=0; l < n_levels; l++){
                help_[l] = local_zeros(n_dofs_[l]);
            }

            return static_cast<Child*>(this)->init_memory_impl(n_dofs_);
            exit(0);


            // // precompute norms of prolongation operators needed for projections of constraints...
            // for(auto l = 0; l < fine_level; l++)
            //     constraints_memory_.P_inf_norm[l] = this->transfer(l).interpolation_inf_norm();

        }



        virtual bool check_feasibility(const SizeType & level, const Vector & x)
        {
            SizeType n_terminates = 0;
            // active lower/upper
            {
                auto d_u = const_device_view(active_upper(level));
                auto d_l = const_device_view(active_lower(level));
                auto d_x = const_device_view(x);


                Device::parallel_reduce(range(x), UTOPIA_LAMBDA(const SizeType i) -> SizeType
                {
                    const Scalar xi = d_x.get(i);
                    const Scalar li = d_l.get(i);
                    const Scalar ui = d_u.get(i);

                    return static_cast<SizeType>(xi < li || xi > ui);
                }, n_terminates);
            }

            bool terminate = n_terminates > 0;
            return x.comm().disjunction(terminate);
        }


      void get_projection(const Vector &lb, const Vector &ub, Vector & x)
      {

        {
            auto d_lb = const_device_view(lb);
            auto d_ub = const_device_view(ub);

            parallel_transform(x, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
            {
                Scalar li = d_lb.get(i);
                Scalar ui = d_ub.get(i);
                if(li >= xi){
                  return li;
                }
                else{
                  return (ui <= xi) ? ui : xi;
                }
            });
        }

      }        


        virtual Scalar criticality_measure_inf(const SizeType & level, const Vector & x, const Vector & g)
        {
            help_[level] = x - g; 
            get_projection(active_lower(level), active_upper(level), help_[level]);
            help_[level] -= x;

            return norm2(help_[level]);
        }


        void init_level(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine)
        {

            // constr_impl_.init_level(level, x_finer_level, x_level, delta_fine); 

            return static_cast<Child*>(this)->init_level_impl(level, x_finer_level, x_level, delta_fine);



            // const SizeType finer_level = level+1;

            // if( IdentityTransfer<Matrix, Vector>* id_transfer =  dynamic_cast<IdentityTransfer<Matrix, Vector>* > (transfer_[level].get()))
            // {
            //     // constraints_memory_.tr_upper[level] =  x_finer_level + local_values(local_size(x_finer_level), delta_fine);
            //     // constraints_memory_.tr_lower[level] =  x_finer_level - local_values(local_size(x_finer_level), delta_fine);

            //     {
            //         auto d_x_finer_level = const_device_view(x_finer_level);

            //         parallel_transform(constraints_memory_.tr_upper[level], UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
            //         {
            //             return d_x_finer_level.get(i) + delta_fine; 
            //         });

            //         parallel_transform(constraints_memory_.tr_lower[level], UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
            //         {
            //             return d_x_finer_level.get(i) - delta_fine; 
            //         });                    
            //     }
            // }
            // else
            // {
            //     //----------------------------------------------------------------------------
            //     //     soft projection of tr bounds
            //     //----------------------------------------------------------------------------
            //     // Vector tr_fine_last_lower = x_finer_level - local_values(local_size(x_finer_level), delta_fine);
            //     // {
            //     //     ReadAndWrite<Vector> rv(tr_fine_last_lower);
            //     //     Read<Vector> rl(constraints_memory_.tr_lower[finer_level]);

            //     //     Range r = range(tr_fine_last_lower);

            //     //     for(SizeType i = r.begin(); i != r.end(); ++i){
            //     //         tr_fine_last_lower.set(i, std::max(constraints_memory_.tr_lower[finer_level].get(i), tr_fine_last_lower.get(i)));
            //     //     }
            //     // }

            //     {
            //         auto d_tr_lb    = const_device_view(constraints_memory_.tr_lower[finer_level]);
            //         auto d_x_finer  = const_device_view(x_finer_level);

            //         parallel_each_write(constraints_memory_.active_lower[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
            //         {
            //             auto val = d_x_finer.get(i) - delta_fine;
            //             auto lbi = d_tr_lb.get(i); 

            //             return std::max(lbi, val);
            //         });
            //     }

            //     // Vector tr_fine_last_upper = x_finer_level + local_values(local_size(x_finer_level), delta_fine);
            //     // {
            //     //     ReadAndWrite<Vector> rv(tr_fine_last_upper);
            //     //     Read<Vector> rl(constraints_memory_.tr_upper[finer_level]);

            //     //     Range r = range(tr_fine_last_upper);

            //     //     for(SizeType i = r.begin(); i != r.end(); ++i){
            //     //         tr_fine_last_upper.set(i, std::min(constraints_memory_.tr_upper[finer_level].get(i), tr_fine_last_upper.get(i)));
            //     //     }
            //     // }

            //     {
            //         auto d_tr_ub    = const_device_view(constraints_memory_.tr_upper[finer_level]);
            //         auto d_x_finer  = const_device_view(x_finer_level);

            //         parallel_each_write(constraints_memory_.active_upper[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
            //         {
            //             auto val = d_x_finer.get(i) + delta_fine;
            //             auto ubi = d_tr_ub.get(i); 

            //             return std::min(ubi, val);
            //         });
            //     }


            //     //------------------------ new version, taking into account  positive and negative elements projection separatelly -----------------
            //     transfer_[level]->project_down_positive_negative(constraints_memory_.active_lower[finer_level], constraints_memory_.active_upper[finer_level], constraints_memory_.tr_lower[level]);
            //     transfer_[level]->project_down_positive_negative(constraints_memory_.active_upper[finer_level], constraints_memory_.active_lower[finer_level], constraints_memory_.tr_upper[level]);

            // }


            // if(has_box_constraints_)
            // {

            //     if( IdentityTransfer<Matrix, Vector>* id_transfer =  dynamic_cast<IdentityTransfer<Matrix, Vector>* > (transfer_[level].get()))
            //     {
            //         // this should be done only on first iteration as it is always same... 
            //         constraints_memory_.x_lower[level] = constraints_memory_.x_lower[finer_level]; 
            //         constraints_memory_.x_upper[level] = constraints_memory_.x_upper[finer_level]; 
            //     }
            //     else
            //     {
            //         //----------------------------------------------------------------------------
            //         //     projection of variable bounds to the coarse level
            //         //----------------------------------------------------------------------------
            //         // Vector lx =  constraints_memory_.x_lower[finer_level] - x_finer_level;
            //         // Scalar lower_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * max(lx);
            //         // constraints_memory_.x_lower[level] = x_level + local_values(local_size(x_level), lower_multiplier);

            //         constraints_memory_.active_lower[finer_level] =  constraints_memory_.x_lower[finer_level] - x_finer_level;
            //         Scalar lower_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * max(constraints_memory_.active_lower[finer_level]);    
            //         {
            //             auto d_x_level    = const_device_view(x_level);
            //             parallel_each_write(constraints_memory_.x_lower[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
            //             {
            //                 auto xi = x_level.get(i);
            //                 return xi + lower_multiplier; 
            //             });
            //         }

            //         // Vector ux =  constraints_memory_.x_upper[finer_level] - x_finer_level;
            //         // Scalar upper_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * min(ux);
            //         // constraints_memory_.x_upper[level] = x_level + local_values(local_size(x_level), upper_multiplier);


            //         constraints_memory_.active_upper[finer_level] =  constraints_memory_.x_lower[finer_level] - x_finer_level;
            //         Scalar upper_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * min(constraints_memory_.active_upper[finer_level]);    
            //         {
            //             auto d_x_level    = const_device_view(x_level);
            //             parallel_each_write(constraints_memory_.x_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
            //             {
            //                 auto xi = x_level.get(i);
            //                 return xi + upper_multiplier; 
            //             });
            //         }


            //     }

            //     //----------------------------------------------------------------------------
            //     //     intersect bounds on the coarse level
            //     //----------------------------------------------------------------------------
            //     // constraints_memory_.active_upper[level] = local_zeros(local_size(x_level));
            //     // constraints_memory_.active_lower[level] = local_zeros(local_size(x_level));
            //     {
            //         // Write<Vector>   rv(constraints_memory_.active_upper[level]), rw(constraints_memory_.active_lower[level]);
            //         // Read<Vector>    rl(x_level), rq(constraints_memory_.x_lower[level]), re(constraints_memory_.x_upper[level]), rr(constraints_memory_.tr_lower[level]), rt(constraints_memory_.tr_upper[level]);

            //         // Range r = range(x_level);

            //         // for(SizeType i = r.begin(); i != r.end(); ++i)
            //         // {
            //         //     constraints_memory_.active_upper[level].set(i, std::min(constraints_memory_.tr_upper[level].get(i), constraints_memory_.x_upper[level].get(i)));
            //         //     constraints_memory_.active_lower[level].set(i, std::max(constraints_memory_.tr_lower[level].get(i), constraints_memory_.x_lower[level].get(i)));
            //         // }

            //         auto d_tr_ub    = const_device_view(constraints_memory_.tr_upper[level]);
            //         auto d_x_ub     = const_device_view(constraints_memory_.x_upper[level]);


            //         parallel_each_write(constraints_memory_.active_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
            //         {
            //             auto val1 = d_tr_ub.get(i);
            //             auto val2 = d_x_ub.get(i); 
            //             return std::min(val1, val2); 
            //         });

            //         auto d_tr_lb    = const_device_view(constraints_memory_.tr_lower[level]);
            //         auto d_x_lb     = const_device_view(constraints_memory_.x_lower[level]);

            //         parallel_each_write(constraints_memory_.active_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
            //         {
            //             auto val1 = d_tr_lb.get(i);
            //             auto val2 = d_x_lb.get(i); 
            //             return std::max(val1, val2); 
            //         });                    
            //     }
            // }
            // else
            // {
            //     constraints_memory_.active_upper[level] = constraints_memory_.tr_upper[level];
            //     constraints_memory_.active_lower[level] = constraints_memory_.tr_lower[level];
            // }


            // base_.interface(); 
            // exit(0);

        }



    protected:
        BoxConstraints box_constraints_;        // constraints on the finest level....
        bool has_box_constraints_;               // as we can run rmtr with inf. norm also without constraints...      

        const std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > & transfer_; 

    protected:
        std::vector<Vector> help_;         

    };






    template<class Matrix, class Vector>
    class IdentityConstraints : public MultilevelVariableBoundSolverInterface<Matrix, Vector, IdentityConstraints<Matrix, Vector> >
    {
        public:
            typedef UTOPIA_SCALAR(Vector)                           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

            typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, IdentityConstraints<Matrix, Vector> > Base; 

            IdentityConstraints(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> & transfer) : Base(transfer)
            {

            }


            void init_memory_impl(const std::vector<SizeType> & n_dofs_)
            {
                std::cout << "IdentityConstraints memory implementation\n";
                constraints_memory_.init_memory(n_dofs_); 
            }

            void init_level_impl(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine)
            {
                std::cout << "IdentityConstraints init_level_impl \n";

                // {
                //     auto d_x_finer_level = const_device_view(x_finer_level);

                //     parallel_transform(constraints_memory_.tr_upper[level], UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
                //     {
                //         return d_x_finer_level.get(i) + delta_fine; 
                //     });

                //     parallel_transform(constraints_memory_.tr_lower[level], UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
                //     {
                //         return d_x_finer_level.get(i) - delta_fine; 
                //     });                    
                // }

                // constraints_memory_.x_lower[level] = constraints_memory_.x_lower[finer_level]; 
                // constraints_memory_.x_upper[level] = constraints_memory_.x_upper[finer_level]; 

                // auto d_tr_ub    = const_device_view(constraints_memory_.tr_upper[level]);
                // auto d_x_ub     = const_device_view(constraints_memory_.x_upper[level]);


                // parallel_each_write(constraints_memory_.active_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                // {
                //     auto val1 = d_tr_ub.get(i);
                //     auto val2 = d_x_ub.get(i); 
                //     return std::min(val1, val2); 
                // });

                // auto d_tr_lb    = const_device_view(constraints_memory_.tr_lower[level]);
                // auto d_x_lb     = const_device_view(constraints_memory_.x_lower[level]);

                // parallel_each_write(constraints_memory_.active_upper[level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                // {
                //     auto val1 = d_tr_lb.get(i);
                //     auto val2 = d_x_lb.get(i); 
                //     return std::max(val1, val2); 
                // });   

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
    };





    template<class Matrix, class Vector>
    class PureTRBounds : public MultilevelVariableBoundSolverInterface<Matrix, Vector, PureTRBounds<Matrix, Vector> >
    {
        public:
            typedef UTOPIA_SCALAR(Vector)                           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

            typedef utopia::MultilevelVariableBoundSolverInterface<Matrix, Vector, PureTRBounds<Matrix, Vector> > Base; 

            PureTRBounds(const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> & transfer) : Base(transfer)
            {

            }


            void init_memory_impl(const std::vector<SizeType> & n_dofs_)
            {
                std::cout << "PureTRBounds memory implementation\n";
                constraints_memory_.init_memory(n_dofs_); 
            }

            void init_level_impl(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine)
            {
                std::cout << "PureTRBounds init_level_impl \n";

                auto finer_level = level + 1; 

                {
                    auto d_x_finer      = const_device_view(x_finer_level);
                    auto d_tr_lb        = const_device_view(constraints_memory_.active_lower[finer_level]);

                    parallel_each_write(this->help_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        auto val = d_x_finer.get(i) - delta_fine;
                        auto lbi = d_tr_lb.get(i); 
                        return std::max(lbi, val);
                    });   
                }

                // only true projection has all elements positive
                this->transfer_[level]->project_down(this->help_[finer_level], constraints_memory_.active_lower[level]);
                

                {
                    auto d_x_finer  = const_device_view(x_finer_level);
                    auto d_tr_ub    = const_device_view(constraints_memory_.active_upper[finer_level]);

                    parallel_each_write(this->help_[finer_level], UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        auto val = d_x_finer.get(i) + delta_fine;
                        auto ubi = d_tr_ub.get(i); 
                        return std::min(ubi, val);
                    });   
                }
                // only true projection has all elements positive
                this->transfer_[level]->project_down(this->help_[finer_level], constraints_memory_.active_upper[level]);                

                // disp(constraints_memory_.active_lower[level], "lower"); 
                // disp(constraints_memory_.active_upper[level], "upper"); 


                //------------------------ we should take into account  positive and negative elements projection separatelly -----------------
                // this->transfer_[level]->project_down_positive_negative(constraints_memory_.active_lower[finer_level], constraints_memory_.active_upper[finer_level], this->help_[level]);
                // this->transfer_[level]->project_down_positive_negative(constraints_memory_.active_upper[finer_level], constraints_memory_.active_lower[finer_level], this->help_[level]);
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
    };




}
#endif //UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP
