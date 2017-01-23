/*
* @Author: alenakopanicakova
* @Date:   2016-05-26
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-07-29
*/

#ifndef UTOPIA_HORIZONTAL_DECOMPOSITION_HPP
#define UTOPIA_HORIZONTAL_DECOMPOSITION_HPP


     namespace utopia 
     {
        /**
         * @brief      Class carries over interpolation/restriction operators/routines for domain decomposition framework. 
         *
         * @tparam     Matrix  
         * @tparam     Vector  
         */
        template <class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
        class HorizontalDecomposition 
        {

        public:

        HorizontalDecomposition()
        {

        }

        ~HorizontalDecomposition(){} 


        virtual bool init(const GlobalVector &x) const = 0; 
        
        virtual bool interpolate(const LocalVector &, GlobalVector &) const = 0;
        virtual bool interpolate(const LocalMatrix &, GlobalMatrix &) const = 0;

        virtual bool restrict(const GlobalVector &, LocalVector &) const = 0;
        virtual bool restrict(const GlobalMatrix &, LocalMatrix &) const = 0;


    protected:        



    };

}

#endif //UTOPIA_HORIZONTAL_DECOMPOSITION_HPP

