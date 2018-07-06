#ifndef UTOPIA_LEVEL_MEMORY_HPP
#define UTOPIA_LEVEL_MEMORY_HPP
#include "utopia_Core.hpp"

namespace utopia 
{

    template<class Matrix, class Vector>
    class FASLevelMemory
    {
        public: 
            void init(const int n_levels)
            {
                x.resize(n_levels);
                x_0.resize(n_levels);

                g.resize(n_levels);                 
                g_diff.resize(n_levels);

                c.resize(n_levels);
            }

        std::vector<Vector> x, x_0, g, g_diff, c; 
    }; 
       

    template<class Matrix, class Vector>
    class RMTRLevelMemory
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;

        public: 
            void init(const int n_levels)
            {
                x.resize(n_levels);
                x_0.resize(n_levels);

                g.resize(n_levels);                 
                g_diff.resize(n_levels);

                s.resize(n_levels);

                H_diff.resize(n_levels);

                delta.resize(n_levels);
            }

        std::vector<Scalar> delta; 
        std::vector<Vector> x, x_0, g, g_diff, s; 
        std::vector<Matrix> H_diff; 
    }; 

    
    template<class Vector>
    class ConstraintsLevelMemory
    {

        public: 
            void init(const int n_levels)
            {
                x_lower.resize(n_levels);
                x_upper.resize(n_levels);
            }

        std::vector<Vector> x_lower, x_upper; 
    }; 



}

#endif //UTOPIA_LEVEL_MEMORY_HPP

