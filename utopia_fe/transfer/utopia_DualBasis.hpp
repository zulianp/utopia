#ifndef UTOPIA_DUAL_BASIS_HPP
#define UTOPIA_DUAL_BASIS_HPP

#include "libmesh/elem.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/fe_base.h"
#include "utopia_libmesh_Utils.hpp"
#include "MortarAssemble.hpp"

namespace utopia {

    //@brief from the paper DUAL QUADRATIC MORTAR FINITE ELEMENT METHODS FOR 3D FINITE DEFORMATION CONTACT∗
    class DualBasis {
    public:

        //transposed trafo
        static bool assemble_local_trafo(
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
                return true;
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
                return true;
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
                return true;
            }

            assert(false && "not implemented");
            return false;
        }

        static void assemble_biorth_weights(
            const libMesh::Elem &el,
            const int el_order,
            libMesh::DenseMatrix<libMesh::Real> &weights)
        {
            const auto dim = el.dim();
            std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(dim, libMesh::Order(el_order));

            const int order = order_for_l2_integral(dim, el, el_order, el, el_order);
            libMesh::QGauss qg(dim, libMesh::Order(order));
            biorth_elem->attach_quadrature_rule(&qg);
            biorth_elem->reinit(&el);
            assemble_biorth_weights(*biorth_elem, weights);
        }

        static void assemble_biorth_weights(const libMesh::FEBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights)
        {
            const auto &test = fe.get_phi();
            const auto &JxW  = fe.get_JxW();

            const uint n_test  = test.size();
            const uint n_qp    = test[0].size();

            libMesh::DenseMatrix<libMesh::Real> elmat;
            elmat.resize(n_test, n_test);
            elmat.zero();

            for(uint i = 0; i < n_test; ++i) {
                for(uint j = 0; j < n_test; ++j) {
                    for(uint qp = 0; qp < n_qp; ++qp) {
                        elmat(i, j) += test[i][qp] * test[j][qp] * JxW[qp];
                    }
                }
            }

            assemble_biorth_weights(elmat, weights);
        }

        static void assemble_biorth_weights(
            const libMesh::Elem &el,
            const int el_order,
            const libMesh::DenseMatrix<libMesh::Real> &trafo,
            libMesh::DenseMatrix<libMesh::Real> &weights)
        {
            const auto dim = el.dim();
            std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(dim, libMesh::Order(el_order));

            const int order = order_for_l2_integral(dim, el, el_order, el, el_order);
            libMesh::QGauss qg(dim, libMesh::Order(order));
            biorth_elem->attach_quadrature_rule(&qg);
            biorth_elem->reinit(&el);
            assemble_biorth_weights(*biorth_elem, trafo, weights);
        }

        static void assemble_biorth_weights(
            const libMesh::FEBase &fe,
            const libMesh::DenseMatrix<libMesh::Real> &trafo,
            libMesh::DenseMatrix<libMesh::Real> &weights)
        {
            const auto &test = fe.get_phi();
            const auto &JxW  = fe.get_JxW();

            const uint n_test  = test.size();
            const uint n_qp    = test[0].size();

            libMesh::DenseMatrix<libMesh::Real> elmat;
            elmat.resize(n_test, n_test);
            elmat.zero();

            std::vector<std::vector<double>> trafo_phi(n_test);;
            for(uint i = 0; i < n_test; ++i) {
                trafo_phi[i].resize(n_qp);

                for(uint qp = 0; qp < n_qp; ++qp) {
                    auto &val = trafo_phi[i][qp];
                    val = 0.;

                    for(uint j = 0; j < n_test; ++j) {
                        val += trafo(i, j) * test[j][qp];
                    }
                }
            }

            for(uint i = 0; i < n_test; ++i) {
                for(uint j = 0; j < n_test; ++j) {
                    for(uint qp = 0; qp < n_qp; ++qp) {
                        elmat(i, j) += trafo_phi[i][qp] * trafo_phi[j][qp] * JxW[qp];
                    }
                }
            }
            
            assemble_biorth_weights(elmat, weights);
        }

        static void assemble_biorth_weights(
            libMesh::DenseMatrix<libMesh::Real> &elmat,
            libMesh::DenseMatrix<libMesh::Real> &weights)
        {
            auto n_test = elmat.n();
            libMesh::DenseVector<libMesh::Real> rhs(n_test);
            rhs.zero();

            libMesh::DenseVector<libMesh::Real> sol(n_test);
            sol.zero();

            libMesh::DenseVector<libMesh::Real> sum_elmat(n_test);
            sum_elmat.zero();

            weights.resize(n_test, n_test);
            weights.zero();

            for(uint i = 0; i < n_test; ++i) {
                for(uint j = 0; j < n_test; ++j) {
                    sum_elmat(i) += elmat(i, j);
                }
            }

            for(uint i = 0; i < n_test; ++i) {
                if(std::abs(sum_elmat(i)) < 1e-16) {
                    sum_elmat(i) = 0;

                    //set identity row where not defined
                    for(uint j = 0; j < n_test; ++j) {
                        elmat(i, j) = (i == j);
                    }
                }
            }

            for(uint i = 0; i < n_test; ++i) {
                if(sum_elmat(i) == 0) {
                    continue;
                }

                rhs(i) = sum_elmat(i);

                elmat.cholesky_solve(rhs, sol);

                for(uint j = 0; j < n_test; ++j) {
                    weights(i, j) = sol(j);
                }

                rhs(i) = 0;
            }

            //normalization for consistently scaled coefficients
            for(uint i = 0; i < n_test; ++i) {
                if(sum_elmat(i) == 0) {
                    continue;
                }

                libMesh::Real t = 0;
                for(uint j = 0; j < n_test; ++j) {
                    t += weights(i, j);
                }

                for(uint j = 0; j < n_test; ++j) {
                    weights(i, j) *= 1./t;
                }
            }
        }
    };
}

#endif //UTOPIA_DUAL_BASIS_HPP

