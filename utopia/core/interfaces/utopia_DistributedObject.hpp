#ifndef UTOPIA_DISTRIBUTED_OBJECT_HPP
#define UTOPIA_DISTRIBUTED_OBJECT_HPP

#include "utopia_Communicator.hpp"

namespace utopia {

	class DistributedObject {
	public:
            virtual ~DistributedObject() = default;
            virtual Communicator &comm() = 0;
            virtual const Communicator &comm() const = 0;
	};

}

#endif //UTOPIA_DISTRIBUTED_OBJECT_HPP
