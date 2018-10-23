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
            using TRSubproblem::tr_constrained_solve;

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
        

            virtual bool tr_constrained_solve(const Operator<Vector> &H, const Vector &g, Vector &s, const BoxConstraints & constraints)
            {
                return false;
            };



    };
}

#endif //TR_BOX_SUBPROBLEM
