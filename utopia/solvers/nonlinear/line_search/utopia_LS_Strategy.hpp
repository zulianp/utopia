/*
* @Author: alenakopanicakova
* @Date:   2016-03-10
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-02
*/ 
#ifndef UTOPIA_LS_STRATEGY_HPP
#define UTOPIA_LS_STRATEGY_HPP
#include "utopia_Parameters.hpp"    
#include "utopia_FunctionNormalEq.hpp"

namespace  utopia 
{   
    /**
     * @brief      Base class for different line-search strategies. 
     *
     * @tparam     Matrix 
     * @tparam     Vector 
     */
    template<class Matrix, class Vector>
    class LSStrategy 
    {

        typedef UTOPIA_SCALAR(Vector) Scalar;
    public:

        LSStrategy(const Parameters  params = Parameters())  
        {
                set_parameters(params); 
        };

        /**
         * @brief      Gets the alpha = step-size. Function needs to be provided by actual LS strategies. 
         *
         * @param      
         * @param[in]  
         *
         * @return     The alpha.
         */
        virtual bool get_alpha(Function<Matrix, Vector> &, const Vector &, const Vector& , const Vector &, Scalar &) = 0;
        virtual bool get_alpha(LeastSquaresFunction<Matrix, Vector> &, const Vector &, const Vector& , const Vector &, Scalar &) = 0;
        virtual bool set_parameters(const Parameters /*params*/){ return true; }
    };
}

#endif //UTOPIA_LS_STRATEGY_HPP
