#ifndef UTOPIA_CONTACT_ASSEMBLER_HPP
#define UTOPIA_CONTACT_ASSEMBLER_HPP

#include "utopia_LibMeshFunctionSpaceAdapter.hpp"
#include "utopia_IContact.hpp"

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

        USparseMatrix B, D, D_inv, Q, Q_inv, T, orthogonal_trafo, complete_transformation;
        UVector weighted_gap, gap;
        UVector weighted_normal, normal;
        UVector area;
        UVector is_contact;
        UVector is_glue_node;
        UVector inv_mass_vector;

        // Factorization<USparseMatrix, UVector> solver;
    };

    class ContactAssembler final : public IContact {
    public:
        using IContact::assemble;

        void read(Input &) override;

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
            const ContactParams &params) override
        {
            return assemble(*mesh, *dof_map, params);
        }

        //retro-compatiblity
        bool init_no_contact(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map) override
        {
            return init_no_contact(*mesh, *dof_map);
        }

        void couple(const UVector &in, UVector &out) const override;
        void uncouple(const UVector &in, UVector &out) const override;
        void couple(const USparseMatrix &in, USparseMatrix &out) const override;

        const UVector &gap() const override;
        UVector &gap() override;

        inline void apply_orthogonal_trafo(const UVector &in, UVector &out) const override
        {
            assert(contact_tensors_);
            out = contact_tensors_->orthogonal_trafo * in;
        }

        inline const USparseMatrix &orthogonal_trafo() const override
        {
            return contact_tensors_->orthogonal_trafo;
        }
            
        inline const UVector &normals() const override { 
            assert(contact_tensors_);
            return contact_tensors_->normal;
        }

        void remove_mass(const UVector &in, UVector &out) const override;

        inline const UVector &is_contact_node() const override { 
            assert(contact_tensors_);
            return contact_tensors_->is_contact;
        }

        inline const UVector &is_glue_node() const override { 
            assert(contact_tensors_);
            return contact_tensors_->is_glue_node;
        }
        
        inline bool initialized() const override
        {
            return static_cast<bool>(contact_tensors_);
        }

        inline bool has_contact() const override
        {
            return has_contact_;
        }

        ContactAssembler() : has_contact_(false) {}

        inline void print_debug_info() override
        {

        }

    private:
        std::shared_ptr<ContactTensors> contact_tensors_;
        bool has_contact_;
    };

}

#endif //UTOPIA_CONTACT_ASSEMBLER_HPP