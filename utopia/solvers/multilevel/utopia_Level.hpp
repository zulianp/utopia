/*
* @Author: alenakopanicakova
* @Date:   2016-03-29
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/

#ifndef UTOPIA_ONE_LEVEL_HPP
#define UTOPIA_ONE_LEVEL_HPP
#include "utopia_Core.hpp"

#include <memory>
#include <assert.h>

     namespace utopia
     {
        /**
         * @brief      Class keeps track on operators.
         *
         *
         * @tparam     Matrix
         * @tparam     Vector
         */
        template<class Matrix, class Vector>
        class Level
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        public:

        Level() : _A(std::make_shared<const Matrix>()) {}
        Level(const std::shared_ptr <const Matrix> & A): _A(A)
        {
            assert(_A);
        }


        virtual ~Level(){}


        /**
         * @brief      Setter of stifness matrix.
         *
         * @param[in]  A     The stifness matrix.
         *
         */
        bool A(const std::shared_ptr <const Matrix> & A)
        {
            assert(A);
            _A = A;
            return true;
        }

        /**
         * @brief      Getter for stifness matrix.
         *
         * @return     The stifness on given level.
         */
        const Matrix &A() const
        {
            assert(_A);
            return *_A;
        }



        /**
         * @brief      Getter for stifness matrix.
         *
         * @return     The stifness on given level.
         */
        std::shared_ptr <const Matrix> A_ptr()
        {
            assert(_A);
            return _A;
        }

        /**
         * @brief      Enforce active set to the system on given level.
         *
         * @param[in]  active_set  The active set (indices).
         * @param      x           The solution vector.
         * @param      b           The right hand side.
         */
        bool enforce_active_set(const std::vector<SizeType> & active_set, Vector & x, Vector & b)
        {
            assert(_A);
            return apply_BC_to_system(_A, x, b, active_set);
        }

    protected:
        std::shared_ptr <const Matrix>  _A;

};

}

#endif //UTOPIA_ONE_LEVEL_HPP

