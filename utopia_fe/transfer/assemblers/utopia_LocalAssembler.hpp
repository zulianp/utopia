#ifndef UTOPIA_LOCAL_ASSEMBLER_HPP
#define UTOPIA_LOCAL_ASSEMBLER_HPP

#include "libmesh/elem.h"
#include "libmesh/fe.h"

namespace utopia {
	class LocalAssembler {
	public:
		using Elem = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;

		virtual ~LocalAssembler() {}

		virtual bool assemble(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			Matrix &mat
			) = 0;
	};
}

#endif //UTOPIA_LOCAL_ASSEMBLER_HPP