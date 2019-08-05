#ifndef UTOPIA_DESCRIBABLE_HPP
#define UTOPIA_DESCRIBABLE_HPP

#include <iostream>

namespace utopia {

    class Describable {
    public:
        virtual ~Describable() {}
        virtual void describe(std::ostream &os = std::cout) const
        {
            os << "Implement me!!!" << std::endl;
        }
    };

}

#endif //UTOPIA_DESCRIBABLE_HPP
