#ifndef UTOPIA_DUAL_BASIS_HPP
#define UTOPIA_DUAL_BASIS_HPP

#include "libmesh/elem.h"
#include "libmesh/dense_matrix.h"

#include "utopia_libmesh_Utils.hpp"

namespace utopia {

    class DualBasis {
    public:

        //assemble for p2
        // static void assemble_local_trafo(
        //     const libMesh::Elem &el,
        //     const int el_order,
        //     const double alpha,
        //     libMesh::DenseMatrix<libMesh::Real> &trafo,
        //     libMesh::DenseMatrix<libMesh::Real> &inv_trafo)
        // {
        //     assert(el_order == 2);
        //     int n = is_tri(el.type()) ? 3 : (is_quad(el.type()) ? 4 : 2);
        //     int n_nodes = el.n_nodes();

        //     assert(n * 2 == n_nodes);

        //     trafo.resize(n_nodes, n_nodes);
        //     inv_trafo.resize(n_nodes, n_nodes);

        //     /////////////////////////////////
        //     for(int i = 0; i < n; ++i) {
        //         //block (0,0)
        //         trafo(i, i)     = 1.0;
        //         inv_trafo(i, i) = 1.0;

        //         //block (1,1)
        //         trafo(n + i, n + i) = (1 - 2*alpha);
        //         inv_trafo(n + i, n + i) = 1./(1 - 2*alpha);

        //     }

        //     /////////////////////////////////
        //     //block (1,0)
        //     for(int i = 0; i < n-1; ++i) {
        //         trafo(n + i, i + 1) = alpha;
        //         inv_trafo(n + i, i + 1) = alpha/(1 - 2.*alpha);
        //     }

        //     trafo(n_nodes-1, n-1) = alpha;
        //     trafo(n_nodes-1, 0)   = alpha;

        //     inv_trafo(n_nodes-1, n-1) = alpha/(1 - 2.*alpha);
        //     inv_trafo(n_nodes-1, 0)   = alpha/(1 - 2.*alpha);
        //     /////////////////////////////////
        // }

        //transposed trafo
        static void assemble_local_trafo(
            const libMesh::ElemType el_type,
            const double alpha,
            libMesh::DenseMatrix<libMesh::Real> &trafo,
            libMesh::DenseMatrix<libMesh::Real> &inv_trafo)
        {
            if(el_type == libMesh::EDGE3) {
                trafo.resize(3, 3);
                inv_trafo.resize(3, 3);

                trafo.zero();
                inv_trafo.zero();

                trafo(0, 0) = 1;
                trafo(1, 1) = 1;
                trafo(0, 2) = alpha;
                trafo(1, 2) = alpha;
                trafo(2, 2) = (1 - 2*alpha);

                /////////////////////////////////////////////

                inv_trafo(0, 0) = 1;
                inv_trafo(1, 1) = 1;
                inv_trafo(0, 2) = alpha/(1 - 2.*alpha);
                inv_trafo(1, 2) = alpha/(1 - 2.*alpha);
                inv_trafo(2, 2) = 1./(1 - 2*alpha);
                return;
            }

            if(el_type == libMesh::TRI6) {
                trafo.resize(6, 6);
                inv_trafo.resize(6, 6);

                trafo(0, 0) = 1;
                trafo(0, 3) = alpha;
                trafo(0, 5) = alpha;

                trafo(1, 1) = 1;
                trafo(1, 3) = alpha;
                trafo(1, 4) = alpha;

                trafo(2, 2) = 1;
                trafo(2, 4) = alpha;
                trafo(2, 5) = alpha;
             
                trafo(3, 3) = (1 - 2*alpha);
                trafo(4, 4) = (1 - 2*alpha);
                trafo(5, 5) = (1 - 2*alpha);

                /////////////////////////////////////////////

                inv_trafo(0, 0) = 1;
                inv_trafo(0, 3) = (1 - 2*alpha);
                inv_trafo(0, 5) = (1 - 2*alpha);

                inv_trafo(1, 1) = 1;
                inv_trafo(1, 3) = (1 - 2*alpha);
                inv_trafo(1, 4) = (1 - 2*alpha);

                inv_trafo(2, 2) = 1;
                inv_trafo(2, 4) = (1 - 2*alpha);
                inv_trafo(2, 5) = (1 - 2*alpha);
                
                inv_trafo(3, 3) = 1./(1 - 2*alpha);
                inv_trafo(4, 4) = 1./(1 - 2*alpha);
                inv_trafo(5, 5) = 1./(1 - 2*alpha);
                return;
            }

            if(el_type == libMesh::QUAD8) {
                trafo.resize(8, 8);
                inv_trafo.resize(8, 8);

                trafo(0, 0) = 1;
                trafo(0, 4) = alpha;
                trafo(0, 7) = alpha;

                trafo(1, 1) = 1;
                trafo(1, 4) = alpha;
                trafo(1, 5) = alpha;

                trafo(2, 2) = 1;
                trafo(2, 5) = alpha;
                trafo(2, 6) = alpha;

                trafo(3, 3) = 1;
                trafo(3, 6) = alpha;
                trafo(3, 7) = alpha;
             
                trafo(4, 4) = (1 - 2*alpha);
                trafo(5, 5) = (1 - 2*alpha);
                trafo(6, 6) = (1 - 2*alpha);
                trafo(7, 7) = (1 - 2*alpha);

                /////////////////////////////////////////////

                inv_trafo(0, 0) = 1;
                inv_trafo(0, 4) = (1 - 2*alpha);
                inv_trafo(0, 7) = (1 - 2*alpha);

                inv_trafo(1, 1) = 1;
                inv_trafo(1, 4) = (1 - 2*alpha);
                inv_trafo(1, 5) = (1 - 2*alpha);

                inv_trafo(2, 2) = 1;
                inv_trafo(2, 5) = (1 - 2*alpha);
                inv_trafo(2, 6) = (1 - 2*alpha);

                inv_trafo(3, 3) = 1;
                inv_trafo(3, 6) = (1 - 2*alpha);
                inv_trafo(3, 7) = (1 - 2*alpha);
                
                inv_trafo(4, 4) = 1./(1 - 2*alpha);
                inv_trafo(5, 5) = 1./(1 - 2*alpha);
                inv_trafo(6, 6) = 1./(1 - 2*alpha);
                inv_trafo(7, 7) = 1./(1 - 2*alpha);
                return;
            }
        }
    };
}

#endif //UTOPIA_DUAL_BASIS_HPP

