#include "utopia_QuadratureUtils.hpp"
#include "MortarAssemble.hpp"

namespace utopia {

    std::shared_ptr<libMesh::QBase> QuadratureUtils::nodal_quad_points(const libMesh::Elem &elem) {
        const int n_nodes = elem.n_nodes();
        const int dim = elem.dim();

        auto q = std::make_shared<QMortar>(dim);
        q->resize(n_nodes);
        auto w = 1. / n_nodes;

        for (int i = 0; i < n_nodes; ++i) {
            q->get_weights()[i] = w;

            switch (dim) {
                case 1: {
                    q->get_points()[i] = libMesh::FE<1, libMesh::LAGRANGE>::inverse_map(&elem, elem.node_ref(i), 1e-10);
                    break;
                }
                case 2: {
                    q->get_points()[i] = libMesh::FE<2, libMesh::LAGRANGE>::inverse_map(&elem, elem.node_ref(i), 1e-10);
                    break;
                }
                case 3: {
                    q->get_points()[i] = libMesh::FE<3, libMesh::LAGRANGE>::inverse_map(&elem, elem.node_ref(i), 1e-10);
                    break;
                }
                default: {
                    assert(false);
                    break;
                }
            }
        }

        return q;
    }

}  // namespace utopia
