#ifndef UTOPIA_FE_I_CONTACT_HPP
#define UTOPIA_FE_I_CONTACT_HPP

#include "utopia.hpp"
#include "utopia_fe_base.hpp"
#include "moonolith_search_radius.hpp"

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
        std::shared_ptr< moonolith::SearchRadius<double> > side_set_search_radius;
        std::vector<std::pair<int, int> > contact_pair_tags;
        std::vector<bool> glued;
        unsigned int variable_number;
        bool use_biorthogonal_basis;

        void describe(std::ostream &os) const;
    };

    class IContact : public Configurable {
    public:
        IContact() {}
        virtual ~IContact() {}

        virtual void read(Input &) override {}

        virtual bool init_no_contact(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map) = 0;

        virtual  bool assemble(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const std::shared_ptr<libMesh::DofMap> &dof_map,
            const ContactParams &params) = 0;

        virtual void couple(const UVector &in, UVector &out) const = 0;
        virtual void couple(const USparseMatrix &in, USparseMatrix &out) const = 0;
        virtual void uncouple(const UVector &in, UVector &out) const = 0;
        
        virtual void apply_orthogonal_trafo(const UVector &in, UVector &out) const = 0;

        virtual const USparseMatrix &orthogonal_trafo() const = 0; 
            
        virtual UVector &gap() = 0;
        virtual const UVector &gap() const  = 0;
        virtual const UVector &normals() const = 0;
        
        virtual void remove_mass(const UVector &in, UVector &out) const = 0;
        virtual const UVector &is_contact_node() const = 0;
        virtual const UVector &is_glue_node() const = 0;

        virtual bool initialized() const = 0;
        virtual bool has_contact() const = 0;
        virtual void print_debug_info() = 0;
    };
}

#endif //UTOPIA_FE_I_CONTACT_HPP

