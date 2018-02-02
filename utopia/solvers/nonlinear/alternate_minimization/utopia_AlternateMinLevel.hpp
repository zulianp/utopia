/*
* @Author: alenakopanicakova
* @Date:   2018-02-02
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2018-02-02
*/

#ifndef UTOPIA_ONE_ALTERNATE_MIN_LEVEL_HPP
#define UTOPIA_ONE_ALTERNATE_MIN_LEVEL_HPP
#include "utopia_Core.hpp"

#include <memory>
#include <assert.h>

    namespace utopia 
    {

        // THIS CLASS NEEDS WAY MUCH BETTER INTERFACE ........... 
        template<class Matrix, class Vector, class FunctionType>
        class AlternateMinLevel 
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
            
            public:

            AlternateMinLevel(const FunctionType & mono, const FunctionType & stag1, const FunctionType & stag2 ) :
            monolithic(mono), 
            stag1(stag1), 
            stag2(stag2)
            {

            }

            virtual ~AlternateMinLevel(){} 






            // protected:        
                std::function< void(const Vector &, Vector &) > from_monolithic_to_stag1; 
                std::function< void(const Vector &, Vector &) > from_monolithic_to_stag2; 
                std::function< void(const Vector &, Vector &) > from_stag1_to_monolithic; 
                std::function< void(const Vector &, Vector &) > from_stag2_to_monolithic; 


                std::function< void(const Vector &) > from_stag1_to_stag2; 
                std::function< void(const Vector &) > from_stag2_to_stag1; 

                FunctionType monolithic; 
                FunctionType stag1; 
                FunctionType stag2; 
        };    

    }

#endif //UTOPIA_ONE_ALTERNATE_MIN_LEVEL_HPP

