#ifndef UTOPIA_MEMORY_BASE_HPP
#define UTOPIA_MEMORY_BASE_HPP

#include "utopia_SolutionStatus.hpp"

namespace utopia
{

    template<class Vector>
    class MemoryInterface {
    public:
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        
        virtual ~MemoryInterface() {}

        virtual void init_memory(const SizeType & ls) = 0; 
    };
}

#endif //UTOPIA_MEMORY_BASE_HPP
