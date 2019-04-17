#ifndef UTOPIA_DUAL_BASIS_HPP
#define UTOPIA_DUAL_BASIS_HPP

namespace utopia {

    class DualBasis {
    public:

        //assemble for p2
        static void assemble_local_trafo(
            const libMesh::Elem &el,
            const int el_order,
            const double alpha,
            libMesh::DenseMatrix<libMesh::Real> &trafo,
            libMesh::DenseMatrix<libMesh::Real> &inv_trafo)
        {
            assert(el_order == 2);
            int n = is_tri(el.type()) ? 3 : (is_quad(el.type()) ? 4 : 2);
            int n_nodes = el.n_nodes();

            assert(n * 2 == n_nodes);

            trafo.resize(n_nodes, n_nodes);
            inv_trafo.resize(n_nodes, n_nodes);

            /////////////////////////////////
            for(int i = 0; i < n; ++i) {
                //block (0,0)
                inv_trafo(i, i) = 1.0;

                //block (1,1)
                inv_trafo(n + i, n + i) = (1 - 2*alpha);
            }

            /////////////////////////////////
            //block (1,0)
            for(int i = 0; i < n-1; ++i) {
                inv_trafo(n + i, i + 1) = alpha;
            }

            inv_trafo(n_nodes-1, n-1) = alpha;
            inv_trafo(n_nodes-1, 0)   = alpha;
            /////////////////////////////////
        }
    };
}

#endif //UTOPIA_DUAL_BASIS_HPP

