#ifndef UTOPIA_CONTACT_ASSEMBLER_HPP
#define UTOPIA_CONTACT_ASSEMBLER_HPP

#include "utopia_LibMeshFunctionSpaceAdapter.hpp"
#include "utopia_Contact.hpp"

#include <cassert>

namespace utopia {

    class ContactTensors {
    public:

        void convert(
            const USparseMatrix &perm,
            const USparseMatrix &vector_perm,
            ContactTensors &out) const;

        
        void finalize(const SizeType spatial_dim, const bool normalize = true);

        static bool check_op(const USparseMatrix &T);

    // private:
        USparseMatrix B_x, D_x, Q_x;

        USparseMatrix B, D, Q, T, orthogonal_trafo, complete_transformation;
        UVector weighted_gap, gap;
        UVector weighted_normal, normal;
        UVector area;
        UVector is_contact;

        Factorization<USparseMatrix, UVector> solver;
    };

    class ContactAssembler {
    public:
        bool assemble(
            libMesh::MeshBase &mesh,
            libMesh::DofMap &dof_map,
            const ContactParams &params);

        void couple(const UVector &in, UVector &out);
        void uncouple(const UVector &in, UVector &out);
        void couple(const USparseMatrix &in, USparseMatrix &out);

        const UVector &gap() const;

    private:
        std::shared_ptr<ContactTensors> contact_tensors_;
    };

}

#endif //UTOPIA_CONTACT_ASSEMBLER_HPP