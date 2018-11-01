#ifndef UTOPIA_QUASI_RMTR_HPP
#define UTOPIA_QUASI_RMTR_HPP

#include "utopia_TRSubproblem.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_RMTR.hpp"


namespace utopia 
{

    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class QuasiRMTR :   public RMTR<Matrix, Vector, CONSISTENCY_LEVEL>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector>     TRSubproblem; 

        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;


        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;


        typedef utopia::BoxConstraints<Vector>                          BoxConstraints;
        typedef utopia::RMTR<Matrix, Vector, CONSISTENCY_LEVEL>         RMTR;

    public:

        QuasiRMTR(  const std::shared_ptr<TRSubproblem> &tr_subproblem_coarse,  
                    const std::shared_ptr<TRSubproblem> &tr_subproblem_smoother,  
                    const Parameters params = Parameters()): 
                    RMTR(tr_subproblem_coarse, tr_subproblem_smoother, params)
        {
            set_parameters(params); 
        }

        virtual ~QuasiRMTR()
        {
            
        } 
        

        void set_parameters(const Parameters params) override
        {
            RMTR::set_parameters(params);    
        }


        virtual std::string name() override 
        { 
            return "Quasi_RMTR";  
        }
        

    protected:
        virtual void init_memory(const SizeType & fine_local_size) override 
        {   
            RMTR::init_memory(fine_local_size);
        }

        virtual void init_coarse_level(const SizeType & level) override
        {
            // init delta on the coarser level...
            this->memory_.delta[level-1]  = this->memory_.delta[level]; 
        }


    protected:   
        ConstraintsLevelMemory <Vector>         constraints_memory_;


    };

}

#endif //UTOPIA_QUASI_RMTR_HPP