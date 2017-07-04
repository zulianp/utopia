/*
* @Author: alenakopanicakova
* @Date:   2017-06-15
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/
#ifndef TR_BOX_SUBPROBLEM
#define TR_BOX_SUBPROBLEM
#include <string>
#include "utopia_TRSubproblem.hpp"
#include "utopia_BoxConstraints.hpp"

namespace  utopia 
{

    /**
     * @brief      Wrapper for TR subproblems 
     */
    template<class Matrix, class Vector>
    class TRBoxSubproblem : public TRSubproblem<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef utopia::TRSubproblem<Matrix, Vector>    TRSubproblem;
        typedef utopia::Preconditioner<Vector>          Preconditioner;
        typedef utopia::BoxConstraints<Vector>          BoxConstraints;

        public:
            TRBoxSubproblem(const Parameters params = Parameters())
                : TRSubproblem(params)
            {
                set_parameters(params); 

            };
            
            virtual ~TRBoxSubproblem( ){}        

            /**
             * @brief      Sets the parameters.
             *
             * @param[in]  params  The parameters
             */
            void set_parameters(const Parameters params) override
            {
                TRSubproblem::set_parameters(params); 
            } 



        public: 
            /**
             * @brief      solve with respect to problem and TR constrians
             *              NEEDS to be provided by sub-classes
             *
             * @param[in]  H             hessian
             * @param[in]  g             graident
             * @param      p_k           new step
             * @param[in]  up_constrain  constrains
             *
             */
            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k, const BoxConstraints & up_constrain) = 0; 


        protected:  


            virtual bool merge_tr_with_pointwise_constrains(const Vector & x_k, const Scalar & radius, const BoxConstraints & pointwise_constrain, Vector & u_f, Vector & l_f)
            {
                if(pointwise_constrain.has_upper_bound())
                {
                    Vector u =  *pointwise_constrain.upper_bound() - x_k; 
                    u_f = local_zeros(local_size(x_k).get(0)); 
                    {   
                        Read<Vector> rv(u); 
                        Write<Vector> wv(u_f); 

                        each_write(u_f, [radius, u](const SizeType i) -> double { 
                            return  (u.get(i) <= radius)  ? u.get(i) : radius; }   );
                    }
                }
                else
                    u_f = radius * local_values(local_size(x_k).get(0), 1.0); ; 


                if(pointwise_constrain.has_lower_bound())
                {
                    Vector l = *(pointwise_constrain.lower_bound()) - x_k; 
                    l_f = local_zeros(local_size(x_k).get(0)); 

                    {   
                        Read<Vector> rv(l); 
                        Write<Vector> wv(l_f); 

                        each_write(l_f, [radius, l](const SizeType i) -> double { 
                            return  (l.get(i) >= -1*radius)  ? l.get(i) : -1 * radius;  }   );
                    }
                }
                else
                    l_f = -1 * radius * local_values(local_size(x_k).get(0), 1.0); ;         

                return true; 
            }



    public: 

        virtual bool  prepare_tr_box_solve(const Scalar & tr_radius, const Vector & x_k, const BoxConstraints & pointwise_constrain, Vector & ub, Vector & lb)
        {
            this->current_radius(tr_radius); 
            merge_tr_with_pointwise_constrains(x_k, tr_radius, pointwise_constrain, ub, lb); 

            return true; 
        }

        
    };
}

#endif //TR_BOX_SUBPROBLEM
