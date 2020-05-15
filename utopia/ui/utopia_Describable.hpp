#ifndef UTOPIA_DESCRIBABLE_HPP
#define UTOPIA_DESCRIBABLE_HPP

#include <iostream>

namespace utopia {

    class Describable {
    public:
        virtual ~Describable() = default;
        virtual void describe(std::ostream &os = std::cout) const { os << "Implement me!!!" << std::endl; }
    };

}  // namespace utopia

#endif  // UTOPIA_DESCRIBABLE_HPP
