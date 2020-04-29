#ifndef UTOPIA_ELEM_HPP
#define UTOPIA_ELEM_HPP

#include <string>

namespace utopia {

    class Elem {
    public:
        Elem() : idx_(-1) {}
        virtual ~Elem() {}

        const std::size_t &idx() const
        {
            return idx_;
        }

        void idx(const std::size_t &idx)
        {
            idx_ = idx;
        }
    private:
        std::size_t idx_;
    };
}

#endif //UTOPIA_ELEM_HPP
