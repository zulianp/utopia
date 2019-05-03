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
        UVector inv_mass_vector;

        // Factorization<USparseMatrix, UVector> solver;
    };

    class ContactAssembler {
    public:
        bool assemble(
            libMesh::MeshBase &mesh,
            libMesh::DofMap &dof_map,
            const ContactParams &params);

        bool init_no_contact(
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map);

        //retro-compatiblity
        inline bool assemble(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map,
            const ContactParams &params)
        {
            return assemble(*mesh, *dof_map, params);
        }

        //retro-compatiblity
        bool init_no_contact(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map)
        {
            return init_no_contact(*mesh, *dof_map);
        }

        void couple(const UVector &in, UVector &out) const;
        void uncouple(const UVector &in, UVector &out) const;
        void couple(const USparseMatrix &in, USparseMatrix &out) const;

        const UVector &gap() const;
        UVector &gap();

        inline void apply_orthogonal_trafo(const UVector &in, UVector &out) const
        {
            assert(contact_tensors_);
            out = contact_tensors_->orthogonal_trafo * in;
        }

        inline const USparseMatrix &orthogonal_trafo() const
        {
            return contact_tensors_->orthogonal_trafo;
        }
        
        
        inline const UVector &normals() const { 
            assert(contact_tensors_);
            return contact_tensors_->normal;
        }
        
        // inline const UVector &inv_mass_vector() const { 
        //     assert(contact_tensors_);
        //     return contact_tensors_->inv_mass_vector;
        // }

        void remove_mass(const UVector &in, UVector &out) const;

        inline const UVector &is_contact_node() const { 
            assert(contact_tensors_);
            return contact_tensors_->is_contact;
        }
        
        inline bool initialized() const
        {
            return static_cast<bool>(contact_tensors_);
        }

        inline bool has_contact() const
        {
            return has_contact_;
        }

        ContactAssembler() : has_contact_(false) {}

    private:
        std::shared_ptr<ContactTensors> contact_tensors_;
        bool has_contact_;
    };

}

#endif //UTOPIA_CONTACT_ASSEMBLER_HPP