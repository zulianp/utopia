#ifndef UTOPIA_FE_CONTACT_HPP
#define UTOPIA_FE_CONTACT_HPP

#include "utopia.hpp"
#include "utopia_fe_base.hpp"

#include <vector>
#include <utility>

namespace libMesh {
    class MeshBase;
    class DofMap;

    namespace Parallel {
        class Communicator;
    }
}

namespace utopia {
    class ContactParams {
    public:
        ContactParams()
        : search_radius(0.1), variable_number(0), use_biorthogonal_basis(true)
        {}

        double search_radius;
        std::vector<std::pair<int, int> > contact_pair_tags;
        unsigned int variable_number;
        bool use_biorthogonal_basis;

        void describe(std::ostream &os) const;
    };

    class Contact {
    public:
        Contact()
        : initialized_(false), has_contact_(false)
        {}

        bool init_no_contact(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map);

        inline bool assemble(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map,
            const ContactParams &params)
        {
            return init(mesh, dof_map, params.search_radius, params.contact_pair_tags, params.variable_number, params.use_biorthogonal_basis);
        }

        void couple(const UVector &in, UVector &out) const;
        void couple(const USparseMatrix &in, USparseMatrix &out) const;
        void uncouple(const UVector &in, UVector &out) const;
        
        inline void apply_orthogonal_trafo(const UVector &in, UVector &out) const
        {
            out = orthogonal_trafo_ * in;
        }

        inline const USparseMatrix &orthogonal_trafo() const
        {
            return orthogonal_trafo_;
        }
            
        inline UVector &gap() { return gap_; }
        inline const UVector &gap() const { return gap_; }

        inline const UVector &normals() const { return normals_; }
        
        // inline const UVector &inv_mass_vector() const { return inv_mass_vector_; }

        inline void remove_mass(const UVector &in, UVector &out) const
        {
            out = e_mul(inv_mass_vector_, in);
        }

        inline const UVector &is_contact_node() const { return is_contact_node_; }
        
        inline bool initialized() const
        {
            return initialized_;
        }

        inline bool has_contact() const
        {
            return has_contact_;
        }

        void print_debug_info();

    private:
        UVector gap_;
        UVector weighted_gap_;
        UVector normals_;
        UVector inv_mass_vector_;
        UVector is_contact_node_;

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

