#ifndef UTOPIA_LOCAL_ASSEMBLER_HPP
#define UTOPIA_LOCAL_ASSEMBLER_HPP

#include "libmesh/elem.h"
#include "libmesh/fe.h"

namespace utopia {



	class LocalAssembler {
	public:
		using Elem   = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;

		enum Type {
			MASTER_X_SLAVE = 0,
			SLAVE_X_SLAVE  = 1
		};

		virtual ~LocalAssembler() {}

		virtual bool assemble(
			const Elem &master,
			FEType master_type,
			const Elem &slave,
			FEType slave_type,
			Matrix &mat
			) = 0;

		virtual bool assemble(
			const Elem &master,
			FEType master_type,
			const Elem &slave,
			FEType slave_type,
			std::vector<Matrix> &mat
			) 
		{
			assert(mat.size() == std::size_t(1));
			return assemble(master, master_type, slave, slave_type, mat[0]);
		}

		virtual int n_forms() const
		{
			return 1;
		}

		virtual Type type(const int index) const
		{
			(void) index;
			return MASTER_X_SLAVE;
		}
	};
}

#endif //UTOPIA_LOCAL_ASSEMBLER_HPP