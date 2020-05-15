#ifndef UTOPIA_OPERATOR_HPP
#define UTOPIA_OPERATOR_HPP

#include "utopia_DistributedObject.hpp"
#include "utopia_Size.hpp"

namespace utopia {

    template <class Vector>
    class Operator : public virtual DistributedObject {
    public:
        using Communicator = typename Traits<Vector>::Communicator;

        ~Operator() override = default;
        virtual bool apply(const Vector &rhs, Vector &sol) const = 0;
        virtual Size size() const = 0;
        virtual Size local_size() const = 0;

        Communicator &comm() override = 0;
        const Communicator &comm() const override = 0;
    };

    template <class Vector, class T>
    const Operator<Vector> &operator_cast(const T &op) {
        return static_cast<const Operator<Vector> &>(op);
    }

}  // namespace utopia

#endif  // UTOPIA_OPERATOR_HPP
