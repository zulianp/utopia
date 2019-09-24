#ifndef UTOPIA_OPERATOR_HPP
#define UTOPIA_OPERATOR_HPP

namespace utopia {

    template<class Vector>
    class Operator {
    public:
        virtual ~Operator() {}
        virtual bool apply(const Vector &rhs, Vector &sol) const = 0;
    };
    
}

#endif //UTOPIA_OPERATOR_HPP
