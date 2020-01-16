#ifndef UTOPIA_BOX_HPP
#define UTOPIA_BOX_HPP

namespace utopia {
    template<class View>
    class Box {
    public:
        View min;
        View max;
    };
}

#endif //UTOPIA_BOX_HPP
