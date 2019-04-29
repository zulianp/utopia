#ifndef UTOPIA_LOCAL_ASSEMBLER_HPP
#define UTOPIA_LOCAL_ASSEMBLER_HPP

#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

#include <cassert>

namespace utopia {

    class LocalContactAssembler {
    public:
        using Elem   = libMesh::Elem;
        using FEType = libMesh::FEType;
        using Matrix = libMesh::DenseMatrix<libMesh::Real>;

        class Result {
        public:
            using Matrix = libMesh::DenseMatrix<libMesh::Real>;
            using Vector = libMesh::DenseVector<libMesh::Real>;

            Matrix coupling_matrix;
            Matrix mass_matrix;
            Matrix normals;
            Vector gap;
        };

        virtual ~LocalContactAssembler() {}

        virtual bool assemble(
            const Elem &master,
            const int master_side,
            FEType master_type,
            const Elem &slave,
            const int slave_side,
            FEType slave_type,
            Result &result
            ) = 0;

        virtual void print_stats(std::ostream &os = std::cout) const
        {
            (void)os;
        }

    };
}

#endif //UTOPIA_LOCAL_ASSEMBLER_HPP
