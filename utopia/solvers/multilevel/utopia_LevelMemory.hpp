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
                s_working.resize(n_levels);

                H_diff.resize(n_levels);

                delta.resize(n_levels);
            }

        std::vector<Scalar> delta; 
        std::vector<Vector> x, x_0, g, g_diff, s, s_working; 
        std::vector<Matrix> H_diff; 
    }; 

    
    template<class Vector>
    class ConstraintsLevelMemory
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;

        public: 
            void init(const int n_levels)
            {
                x_lower.resize(n_levels);
                x_upper.resize(n_levels);

                tr_lower.resize(n_levels);
                tr_upper.resize(n_levels);
                
                active_lower.resize(n_levels);
                active_upper.resize(n_levels);                                

                P_inf_norm.resize(n_levels); 
            }

        std::vector<Vector> x_lower, x_upper, tr_lower, tr_upper, active_lower, active_upper; 
        std::vector<Scalar> P_inf_norm; 
    }; 



}

#endif //UTOPIA_LEVEL_MEMORY_HPP

