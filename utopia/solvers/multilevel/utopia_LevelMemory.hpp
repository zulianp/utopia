#ifndef UTOPIA_LEVEL_MEMORY_HPP
#define UTOPIA_LEVEL_MEMORY_HPP
#include "utopia_Core.hpp"

namespace utopia
{

    enum MultiLevelCoherence{   FIRST_ORDER     = 1,
                                FIRST_ORDER_DF  = 3,
                                SECOND_ORDER    = 2,
                                GALERKIN        = 0};

    template<MultiLevelCoherence T, MultiLevelCoherence U>
    struct is_same : std::false_type {};

    template<MultiLevelCoherence T>
    struct is_same<T, T> : std::true_type {};                                


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
        typedef UTOPIA_SCALAR(Vector)       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;

        public:
            void init_memory(const std::vector<SizeType> & n_dofs_)
            {
                const auto n_levels = n_dofs_.size(); 

                x.resize(n_levels);
                x_0.resize(n_levels);

                s.resize(n_levels);
                s_working.resize(n_levels);
                help.resize(n_levels); 

                delta.resize(n_levels);
                energy.resize(n_levels); 
                gnorm.resize(n_levels); 

                for(auto l=0; l < n_levels; l++){
                    x[l]            = local_zeros(n_dofs_[l]); 
                    x_0[l]          = local_zeros(n_dofs_[l]); 
                    s[l]            = local_zeros(n_dofs_[l]); 
                    s_working[l]    = local_zeros(n_dofs_[l]); 
                    help[l]         = local_zeros(n_dofs_[l]); 
                }
            }            

        std::vector<Scalar> delta, energy, gnorm;
        std::vector<Vector> x, x_0, s, s_working, help;
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

