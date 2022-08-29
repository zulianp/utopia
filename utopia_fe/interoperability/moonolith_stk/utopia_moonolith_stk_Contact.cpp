#include "utopia_moonolith_stk_Contact.hpp"

#include "moonolith_contact.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_moonolith_stk.hpp"

#include "moonolith_affine_transform.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

namespace utopia {

    namespace stk {

        class Contact::Impl : public Configurable {
        public:
            using Mesh_t = utopia::moonolith::Mesh;
            using FunctionSpace_t = utopia::moonolith::FunctionSpace;
            using Contact_t = utopia::moonolith::Contact;

            bool assemble(const FunctionSpace &in_space) {
                // FIXME handle trace spaces
                // if (in_space.mesh().spatial_dimension() > in_space.mesh().manifold_dimension()) {
                //     convert_function_space(in_space, space);
                // } else {
                extract_trace_space(in_space, space);
                // }

                // FIXME this is a hack, we should do this with stk by defining parts
                using IndexArray = typename Traits<FunctionSpace>::IndexArray;
                auto banned_nodes = std::make_shared<IndexArray>();
                in_space.create_boundary_node_list(*banned_nodes);
                contact.set_banned_nodes(banned_nodes);

                return contact.assemble(space);
            }

            void read(Input &in) override {
                contact.read(in);
                in.get("remove_constrained_dofs", remove_constrained_dofs);
            }

            FunctionSpace_t space;
            Contact_t contact;
            bool remove_constrained_dofs{false};
        };

        Contact::Contact() : impl_(utopia::make_unique<Impl>()) {}
        Contact::~Contact() {}

        void Contact::read(Input &in) { impl_->contact.read(in); }
        void Contact::describe(std::ostream &os) const { impl_->contact.describe(os); }

        bool Contact::assemble(FunctionSpace &space) {
            UTOPIA_TRACE_REGION_BEGIN("Contact::assemble");

            bool ok = impl_->assemble(space);

            if (ok) {
                if (impl_->remove_constrained_dofs) {
                    space.apply_constraints(*impl_->contact.complete_transformation(), 0);
                }
            }

            UTOPIA_TRACE_REGION_END("Contact::assemble");
            return ok;
        }

        const Contact::Vector &Contact::gap() const { return impl_->contact.gap(); }
        const Contact::Vector &Contact::is_contact() const { return impl_->contact.is_contact(); }
        const Contact::Vector &Contact::normals() const { return impl_->contact.normals(); }

        void Contact::set_params(const Params &params) { impl_->contact.set_params(params); }

        void Contact::transform(const Matrix &in, Matrix &out) { impl_->contact.transform(in, out); }

        void Contact::transform(const Vector &in, Vector &out) { impl_->contact.transform(in, out); }

        void Contact::inverse_transform(const Vector &in, Vector &out) { impl_->contact.inverse_transform(in, out); }

        std::shared_ptr<Contact::Matrix> Contact::orthogonal_transformation() {
            return impl_->contact.orthogonal_transformation();
        }

        std::shared_ptr<Contact::Matrix> Contact::mass_matrix() { return impl_->contact.mass_matrix(); }

    }  // namespace stk
}  // namespace utopia
