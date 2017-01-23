/*
* @Author: alenakopanicakova
* @Date:   2016-05-18
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2016-09-05
*/

#ifndef UTOPIA_GLOBALANDLOCALFUNCTION_DD_HPP
#define UTOPIA_GLOBALANDLOCALFUNCTION_DD_HPP
#include "utopia_Base.hpp"

namespace utopia 
{

    template <class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
    /**
     * @brief      This function contains informations on global and local nonlinear evaluation routines. 
     *             At the moment, also informations on interpolation and restriction are here. 
     *             In future, this is going to be deleted from here. Please, do not use it anymore. 
     */
    class GLFunction : public Function<GlobalMatrix, GlobalVector> 
    {
            DEF_UTOPIA_SCALAR(GlobalMatrix)

        public:
            virtual ~GLFunction(){}

            virtual bool init(const GlobalVector &x)                        const = 0; 

            // global evaluation things
            virtual bool value(const GlobalVector &, Scalar &)              const = 0;
            virtual bool gradient(const GlobalVector &, GlobalVector &)     const = 0;
            virtual bool hessian(const GlobalVector &, GlobalMatrix &)      const = 0;

            // local evaluation things
            virtual bool local_value(const LocalVector &, Scalar &)         const = 0;
            virtual bool local_gradient(const LocalVector &, LocalVector &) const = 0;
            virtual bool local_hessian(const LocalVector &, LocalMatrix &)  const = 0;

            // decomposition staff
            virtual bool interpolate(const LocalVector &, GlobalVector &)   const = 0;
            virtual bool interpolate(const LocalMatrix &, GlobalMatrix &)   const = 0;

            virtual bool restrict(const GlobalVector &, LocalVector &)      const = 0;
            virtual bool restrict(const GlobalMatrix &, LocalMatrix &)      const = 0;


    };
}
#endif //UTOPIA_GLOBALANDLOCALFUNCTION_DD_HPP


