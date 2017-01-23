/*
* @Author: alenakopanicakova
* @Date:   2016-07-06
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-07-29
*/

#ifdef WITH_PETSC
#ifndef UTOPIA_UTOPIA_LOCAL_REDISTRIBUTE_HPP
#define UTOPIA_UTOPIA_LOCAL_REDISTRIBUTE_HPP

#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"
#include <string>


namespace utopia 
{   
    /**
     * @brief  
     *
     * @tparam     Left   The vector with stored values.
     * @tparam     Right  The vector defining shape. 
     * @tparam     Result  The left vector with local distribution of right vector. 
     */
    template<class Left, class Right>
    class LocalRedistribute : public Expression< LocalRedistribute<Left, Right> > 
    {
    public:
        LocalRedistribute(const Left &left, const Right &right)
                : _left(left), _right(right)
        {}

        inline const Left &left() const 
        {
            return _left;
        }

        inline const Right &right() const 
        {
            return _right;
        }

        std::string getClass() const override
        {
            return "LocalRedistribute<" + left().getClass() + ", " + right().getClass() + ">";
        }


    private:
        UTOPIA_STORE_CONST(Left)  _left;
        UTOPIA_STORE_CONST(Right) _right;
    };

    // template<class Left, class Right>
    // inline LocalRedistribute<Left, Right> local_redistribute(const Expression<Left> &left, const Expression<Right> &right)
    // {
    //     return LocalRedistribute<Left, Right>(left.derived(), right.derived());
    // }

    /**
     * @ingroup  parallel_expressions
     * @brief      Reshapes local sizes of vector x, according to the distribution of shape vector, while preserving data of the original vector.

     * @note       Global sizes of vectors x and shape vector have to be same. 
     *
     * @param[in]  left          The vector to be reshaped. 
     * @param[in]  shape_vector  The vector defining new shape. 
     *
     * @tparam     Left          
     * @tparam     Right         
     *
     * @return     
     */
    template<class Left, class Right>
    inline LocalRedistribute<Left, Right> redistribute_as(const Expression<Left> &left, const Expression<Right> &shape_vector)
    {
        return LocalRedistribute<Left, Right>(left.derived(), shape_vector.derived());
    }


}

#endif //UTOPIA_UTOPIA_LOCAL_REDISTRIBUTE_HPP
#endif //WITH_PETSC
