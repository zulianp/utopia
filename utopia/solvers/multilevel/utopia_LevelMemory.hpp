#ifndef UTOPIA_LEVEL_MEMORY_HPP
#define UTOPIA_LEVEL_MEMORY_HPP
#include "utopia_Core.hpp"

namespace utopia
{

    enum MultiLevelCoherence{   FIRST_ORDER         = 1,
                                FIRST_ORDER_DF      = 3,
                                FIRST_ORDER_MGOPT   = 4,
                                
                                SECOND_ORDER    = 2,
                                GALERKIN        = 0};

    template<MultiLevelCoherence T, MultiLevelCoherence U>
    struct is_same : std::false_type {};

    template<MultiLevelCoherence T>
    struct is_same<T, T> : std::true_type {};         
    
    template<MultiLevelCoherence T, MultiLevelCoherence... Rest>
    struct is_any : std::false_type {};

    template<MultiLevelCoherence T, MultiLevelCoherence First>
    struct is_any<T, First> : is_same<T, First> {};          

    template<MultiLevelCoherence T, MultiLevelCoherence First, MultiLevelCoherence... Rest>
    struct is_any<T, First, Rest...> : std::integral_constant<bool, is_same<T, First>::value || is_any<T, Rest...>::value>
    {};

    // helper function, should be implemented in c++ 14
    template< bool B, class T = void >
    using enable_if_t = typename std::enable_if<B,T>::type;



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
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

        public:
            void init_memory(const std::vector<SizeType> & n_dofs_)
            {
                const auto n_levels = n_dofs_.size(); 

                x_lower.resize(n_levels);
                x_upper.resize(n_levels);

                tr_lower.resize(n_levels);
                tr_upper.resize(n_levels);

                active_lower.resize(n_levels);
                active_upper.resize(n_levels);

                help.resize(n_levels); 

                for(auto l=0; l < n_levels; l++){
                    x_lower[l]  = local_zeros(n_dofs_[l]); 
                    x_upper[l]  = local_zeros(n_dofs_[l]); 

                    tr_lower[l]     = local_zeros(n_dofs_[l]); 
                    tr_upper[l]     = local_zeros(n_dofs_[l]); 

                    active_lower[l] = local_zeros(n_dofs_[l]); 
                    active_upper[l] = local_zeros(n_dofs_[l]); 

                    help[l] = local_zeros(n_dofs_[l]); 
                }
            }

        std::vector<Vector> x_lower, x_upper, tr_lower, tr_upper, active_lower, active_upper, help;
    };


    template<class Vector>
    class ActiveConstraintsLevelMemory
    {
        typedef UTOPIA_SCALAR(Vector)        Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)     SizeType;

        public:
            void init_memory(const std::vector<SizeType> & n_dofs_)
            {
                const Scalar inf = std::numeric_limits<Scalar>::infinity();
                const auto n_levels = n_dofs_.size(); 
                active_lower.resize(n_levels);
                active_upper.resize(n_levels);

                for(auto l=0; l < n_levels; l++){
                    active_lower[l] = local_values(n_dofs_[l], -inf); 
                    active_upper[l] = local_values(n_dofs_[l], inf); 
                }
            }
            
        std::vector<Vector> active_lower, active_upper; 
    };


}

#endif //UTOPIA_LEVEL_MEMORY_HPP



