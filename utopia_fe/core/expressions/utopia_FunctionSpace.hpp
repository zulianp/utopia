#ifndef UTOPIA_FUNCTION_SPACE_HPP
#define UTOPIA_FUNCTION_SPACE_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include <string>

namespace utopia {

    template<class Derived>
    class FunctionSpace {
    public:
        FunctionSpace()
        : subspace_id_(0)
        {}

        inline int subspace_id() const
        {
            return subspace_id_;
        }

        void set_subspace_id(const int id)
        {
            subspace_id_ = id;
        }

        DERIVED_CRT(Derived);
        CONST_DERIVED_CRT(Derived);

        virtual std::string get_class() const
        {
            return "FunctionSpace";
        }

    private:
        int subspace_id_;
    };

    template<class Derived>
    class Traits< FunctionSpace<Derived> > : public Traits<Derived> {};
}

#endif //UTOPIA_FUNCTION_SPACE_HPP
