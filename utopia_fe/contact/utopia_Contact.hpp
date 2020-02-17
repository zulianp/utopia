#ifndef UTOPIA_FE_CONTACT_HPP
#define UTOPIA_FE_CONTACT_HPP

#include "utopia.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_IContact.hpp"

#include <vector>
#include <utility>

namespace utopia {

    class Contact final : public IContact {
    public:

        using IContact::assemble;

        Contact()
        : initialized_(false), has_contact_(false)
        {}

        bool init_no_contact(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map) override;

        inline bool assemble(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map,
            const ContactParams &params) override
        {
            return init(mesh, dof_map, params.search_radius, params.contact_pair_tags, params.variable_number, params.use_biorthogonal_basis);
        }

        void couple(const UVector &in, UVector &out) const override;
        void couple(const USparseMatrix &in, USparseMatrix &out) const override;
        void uncouple(const UVector &in, UVector &out) const override;
        
        inline void apply_orthogonal_trafo(const UVector &in, UVector &out) const override
        {
            out = orthogonal_trafo_ * in;
        }

        inline const USparseMatrix &orthogonal_trafo() const override
        {
            return orthogonal_trafo_;
        }
            
        inline UVector &gap() override { return gap_; }
        inline const UVector &gap() const override { return gap_; }

        inline const UVector &normals() const override { return normals_; }

        inline void remove_mass(const UVector &in, UVector &out) const override
        {
            out = e_mul(inv_mass_vector_, in);
        }

        inline const UVector &is_contact_node() const override { return is_contact_node_; }
        inline const UVector &is_glue_node() const override { return is_glue_node_; }
        
        inline bool initialized() const override
        {
            return initialized_;
        }

        inline bool has_contact() const override
        {
            return has_contact_;
        }

        void print_debug_info() override; 

    private:
        UVector gap_;
        UVector weighted_gap_;
        UVector normals_;
        UVector inv_mass_vector_;
        UVector is_contact_node_;
        UVector is_glue_node_;

        USparseMatrix coupling_;
        USparseMatrix inv_mass_matrix_;
        USparseMatrix transfer_operator_;
        USparseMatrix orthogonal_trafo_;

        USparseMatrix complete_transformation_;
        bool initialized_;
        bool has_contact_;


        bool init(const std::shared_ptr<libMesh::MeshBase> &mesh,
                  const std::shared_ptr<libMesh::DofMap> &dof_map,
                  const double search_radius,
                  const std::vector<std::pair<int, int> > &contact_pair_tags,
                  unsigned int variable_number = 0,
                  const bool use_biorthogonal_basis = true);
    };
}

#endif //UTOPIA_FE_CONTACT_HPP

