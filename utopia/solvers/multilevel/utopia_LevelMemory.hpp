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
        
        


}

#endif //UTOPIA_LEVEL_MEMORY_HPP

