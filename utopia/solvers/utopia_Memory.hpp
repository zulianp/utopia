#ifndef UTOPIA_MEMORY_BASE_HPP
#define UTOPIA_MEMORY_BASE_HPP

#include "utopia_SolutionStatus.hpp"
#include "utopia_Layout.hpp"

namespace utopia
{

    template<class Vector>
    class MemoryInterface {
    public:
        using SizeType     = typename Traits<Vector>::SizeType;
        using Communicator = typename Traits<Vector>::Communicator;
        using Layout       = typename Traits<Vector>::Layout;

        virtual ~MemoryInterface() {}

        virtual void init_memory(const Layout &l) = 0;
    };
}

#endif //UTOPIA_MEMORY_BASE_HPP
