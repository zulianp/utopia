#ifndef UTOPIA_DISTRIBUTED_OBJECT_HPP
#define UTOPIA_DISTRIBUTED_OBJECT_HPP

#include "utopia_Communicator.hpp"

namespace utopia {

    class DistributedObject {
    public:
        virtual ~DistributedObject() = default;
        virtual Communicator &comm() {
            return const_cast<Communicator &>(const_cast<DistributedObject *>(this)->comm());
        }

        virtual const Communicator &comm() const = 0;
    };

}  // namespace utopia

#endif  // UTOPIA_DISTRIBUTED_OBJECT_HPP
