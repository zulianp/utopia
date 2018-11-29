#ifndef UTOPIA_L2_LOCAL_ASSEMBLER_HPP
#define UTOPIA_L2_LOCAL_ASSEMBLER_HPP

#include "utopia_LocalContactAssembler.hpp"

#include "libmesh/elem.h"
#include "libmesh/fe.h"

#include <cassert>
#include <ostream>

namespace utopia {

	class L2LocalContactAssembler final : public LocalContactAssembler {
	public:
		using Elem   = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;

		~L2LocalContactAssembler() {}

		bool assemble(
			const Elem &master,
			const int master_side,
			FEType master_type,
			const Elem &slave,
			const int slave_side,
			FEType slave_type,
			Result &result
			) override;

		void print_stats(std::ostream &os = std::cout) const override;
		
	};
}

#endif //UTOPIA_L2_LOCAL_ASSEMBLER_HPP
