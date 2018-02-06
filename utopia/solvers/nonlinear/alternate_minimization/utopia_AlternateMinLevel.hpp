/*
* @Author: alenakopanicakova
* @Date:   2018-02-02
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2018-02-06
*/

#ifndef UTOPIA_ONE_ALTERNATE_MIN_LEVEL_HPP
#define UTOPIA_ONE_ALTERNATE_MIN_LEVEL_HPP
#include "utopia_Core.hpp"

#include <memory>
#include <assert.h>

    namespace utopia 
    {

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



            AlternateMinLevel(  const FunctionType & mono, 
                                const FunctionType & stag1, 
                                const FunctionType & stag2, 
                                const std::function< void(const Vector &, Vector &) > from_monolithic_to_stag1,
                                const std::function< void(const Vector &, Vector &) > from_monolithic_to_stag2, 
                                const std::function< void(const Vector &, Vector &) > from_stag1_to_monolithic, 
                                const std::function< void(const Vector &, Vector &) > from_stag2_to_monolithic, 
                                const std::function< void(const Vector &) > from_stag1_to_stag2, 
                                const std::function< void(const Vector &) > from_stag2_to_stag1 ):
            monolithic(mono), 
            stag1(stag1), 
            stag2(stag2), 
            from_monolithic_to_stag1(from_monolithic_to_stag1), 
            from_monolithic_to_stag2(from_monolithic_to_stag2), 
            from_stag1_to_monolithic(from_stag1_to_monolithic), 
            from_stag2_to_monolithic(from_stag2_to_monolithic), 
            from_stag1_to_stag2(from_stag1_to_stag2), 
            from_stag2_to_stag1(from_stag2_to_stag1)
            {

            }

            virtual ~AlternateMinLevel(){} 



            virtual FunctionType & fun_monolithic()
            {
                return monolithic; 
            }


            virtual FunctionType & fun_stag1()
            {
                return stag1; 
            }

            virtual FunctionType & fun_stag2()
            {
                return stag2; 
            }


            virtual void transfer_from_monolithic_to_stag1(const Vector & from, Vector & to)
            {
                from_monolithic_to_stag1(from, to); 
            }

            virtual void transfer_from_monolithic_to_stag2(const Vector & from, Vector & to)
            {
                from_monolithic_to_stag2(from, to); 
            }

            virtual void transfer_from_stag1_to_monolithic(const Vector & from, Vector & to)
            {
                from_stag1_to_monolithic(from, to); 
            }

            virtual void transfer_from_stag2_to_monolithic(const Vector & from, Vector & to)
            {
                from_stag2_to_monolithic(from, to); 
            }


            virtual void transfer_from_stag1_to_stag2(const Vector& from)
            {
                from_stag1_to_stag2(from); 
            }

            virtual void transfer_from_stag2_to_stag1(const Vector& from)
            {
                from_stag2_to_stag1(from); 
            }



            virtual void set_transfers( const std::function< void(const Vector &, Vector &) > from_monolithic_to_stag1_in,
                                const std::function< void(const Vector &, Vector &) > from_monolithic_to_stag2_in, 
                                const std::function< void(const Vector &, Vector &) > from_stag1_to_monolithic_in, 
                                const std::function< void(const Vector &, Vector &) > from_stag2_to_monolithic_in, 
                                const std::function< void(const Vector &) > from_stag1_to_stag2_in, 
                                const std::function< void(const Vector &) > from_stag2_to_stag1_in)
            {
                from_monolithic_to_stag1 = from_monolithic_to_stag1_in; 
                from_monolithic_to_stag2 = from_monolithic_to_stag2_in; 
                from_stag1_to_monolithic = from_stag1_to_monolithic_in; 
                from_stag2_to_monolithic = from_stag2_to_monolithic_in; 
                from_stag1_to_stag2      = from_stag1_to_stag2_in; 
                from_stag2_to_stag1      = from_stag2_to_stag1_in; 
            }


            virtual void set_transfer_from_monolithic_to_stag1(const std::function< void(const Vector &, Vector &) > from_monolithic_to_stag1_in )
            {
                from_monolithic_to_stag1 = from_monolithic_to_stag1_in; 
            }

            virtual void set_transfer_from_monolithic_to_stag2(const std::function< void(const Vector &, Vector &) > from_monolithic_to_stag2_in )
            {
                from_monolithic_to_stag2 = from_monolithic_to_stag2_in; 
            }

            virtual void set_transfer_from_stag1_to_monolithic(const std::function< void(const Vector &, Vector &) > from_stag1_to_monolithic_in )
            {
                from_stag1_to_monolithic = from_stag1_to_monolithic_in; 
            }

            virtual void set_transfer_from_stag2_to_monolithic(const std::function< void(const Vector &, Vector &) > from_stag2_to_monolithic_in )
            {
                from_stag2_to_monolithic = from_stag2_to_monolithic_in; 
            }                        

            virtual void set_transfer_from_stag1_to_stag2(const std::function< void(const Vector &) > from_stag1_to_stag2_in )
            {
                from_stag1_to_stag2 = from_stag1_to_stag2_in; 
            }                        

            virtual void set_transfer_from_stag2_to_stag1(const std::function< void(const Vector &) > from_stag2_to_stag1_in )
            {
                from_stag2_to_stag1 = from_stag2_to_stag1_in; 
            }                                                

            protected:        
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

