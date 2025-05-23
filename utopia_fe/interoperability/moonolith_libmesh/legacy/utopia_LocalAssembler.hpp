#ifndef UTOPIA_LOCAL_ASSEMBLER_HPP
#define UTOPIA_LOCAL_ASSEMBLER_HPP

#include "libmesh/elem.h"
#include "libmesh/fe.h"

#include <cassert>

namespace utopia {

    class LocalAssembler {
    public:
        using Elem = libMesh::Elem;
        using FEType = libMesh::FEType;
        using Matrix = libMesh::DenseMatrix<libMesh::Real>;

        enum Type { MASTER_X_SLAVE = 0, SLAVE_X_SLAVE = 1, MASTER_X_MASTER = 3, SLAVE_X_MASTER = 4 };

        virtual ~LocalAssembler() {}

        virtual bool assemble(const Elem &master,
                              FEType master_type,
                              const Elem &slave,
                              FEType slave_type,
                              Matrix &mat) = 0;

        virtual bool assemble(const Elem &master,
                              FEType master_type,
                              const Elem &slave,
                              FEType slave_type,
                              std::vector<Matrix> &mat) {
            assert(n_forms() == int(1));

            if (mat.empty()) {
                mat.resize(1);
            }

            return assemble(master, master_type, slave, slave_type, mat[0]);
        }

        virtual int n_forms() const { return 1; }

        virtual Type type(const int index) const {
            (void)index;
            return MASTER_X_SLAVE;
        }

        virtual void print_stats(std::ostream &os = std::cout) const { (void)os; }
    };
}  // namespace utopia

#endif  // UTOPIA_LOCAL_ASSEMBLER_HPP
