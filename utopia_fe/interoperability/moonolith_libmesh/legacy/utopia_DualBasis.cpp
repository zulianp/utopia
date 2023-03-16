#include "utopia_DualBasis_impl.hpp"

#include "MortarAssemble.hpp"
#include "libmesh/dense_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/reference_elem.h"
#include "utopia_ElementDofMap.hpp"

namespace utopia {

    bool is_diag(const libMesh::DenseMatrix<double> &d, const bool verbose) {
        bool ret = true;
        for (int i = 0; i < d.m(); ++i) {
            for (int j = 0; j < d.n(); ++j) {
                if (i != j) {
                    if (!approxeq(d(i, j), 0.0, 1e-10)) {
                        ret = false;
                        break;
                    }
                }
            }

            if (!ret) break;
        }

        if (verbose && !ret) {
            std::cerr << "---------------\n";
            d.print();
            std::cerr << "---------------\n";
        }

        return ret;
    }

    //@brief from the paper DUAL QUADRATIC MORTAR FINITE ELEMENT METHODS FOR 3D FINITE DEFORMATION CONTACT∗

    void DualBasis::init(const libMesh::ElemType &type, const double alpha) {
        build_trafo_and_weights(type, order, alpha, trafo_, inv_trafo_, weights_);
    }

    void DualBasis::compute_values(const libMesh::FEBase &fe) {
        if (compute_phi) {
            const auto &phi = fe.get_phi();

            uint n_fun = phi.size();
            uint n_qp = phi[0].size();

            phi_.resize(n_fun);

            for (uint p = 0; p < n_fun; ++p) {
                phi_[p].resize(n_qp);

                for (uint k = 0; k < n_qp; ++k) {
                    auto f = 0.0;

                    for (uint i = 0; i < n_fun; ++i) {
                        f += weights_(p, i) * phi[i][k];
                    }

                    phi_[p][k] = f;
                }
            }
        }

        if (compute_dphi) {
            const auto &dphi = fe.get_dphi();

            uint n_fun = dphi.size();
            uint n_qp = dphi[0].size();

            dphi_.resize(n_fun);

            for (uint p = 0; p < n_fun; ++p) {
                dphi_[p].resize(n_qp);

                for (uint k = 0; k < n_qp; ++k) {
                    auto f = weights_(p, 0) * dphi[0][k];

                    for (uint i = 1; i < n_fun; ++i) {
                        f += weights_(p, i) * dphi[i][k];
                    }

                    dphi_[p][k] = f;
                }
            }
        }
    }

    bool DualBasis::build_trafo_and_weights(const libMesh::ElemType type,
                                            const int order,
                                            const double alpha,
                                            libMesh::DenseMatrix<libMesh::Real> &trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &inv_trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &weights) {
        const bool is_higher_order_tp = type == libMesh::HEX27 || type == libMesh::QUAD9;

        if (order != 1 && !is_higher_order_tp) {
            if (!DualBasis::assemble_local_trafo(type, alpha, trafo, inv_trafo)) {
                assert(false);
                return false;
            }

            const auto &ref_elem = libMesh::ReferenceElem::get(type);

            DualBasis::assemble_biorth_weights(ref_elem, order, trafo, weights, false);

            weights.right_multiply(trafo);
        } else {
            const auto &ref_elem = libMesh::ReferenceElem::get(type);

            assemble_biorth_weights(ref_elem, order, weights, !is_higher_order_tp);

            auto n = weights.n();
            trafo.resize(n, n);
            inv_trafo.resize(n, n);

            for (uint i = 0; i < n; ++i) {
                trafo(i, i) = 1.0;
                inv_trafo(i, i) = 1.0;
            }
        }
        return true;
    }

    // transposed trafo
    bool DualBasis::assemble_local_trafo(const libMesh::ElemType el_type,
                                         const double alpha,
                                         libMesh::DenseMatrix<libMesh::Real> &trafo,
                                         libMesh::DenseMatrix<libMesh::Real> &inv_trafo) {
        if (el_type == libMesh::EDGE3) {
            trafo.resize(3, 3);
            inv_trafo.resize(3, 3);

            trafo.zero();
            inv_trafo.zero();

            trafo(0, 0) = 1;
            trafo(1, 1) = 1;
            trafo(0, 2) = alpha;
            trafo(1, 2) = alpha;
            trafo(2, 2) = (1 - 2 * alpha);

            /////////////////////////////////////////////

            inv_trafo(0, 0) = 1;
            inv_trafo(1, 1) = 1;
            inv_trafo(0, 2) = alpha / (1 - 2. * alpha);
            inv_trafo(1, 2) = alpha / (1 - 2. * alpha);
            inv_trafo(2, 2) = 1. / (1 - 2 * alpha);
            return true;
        }

        if (el_type == libMesh::TRI6) {
            trafo.resize(6, 6);
            inv_trafo.resize(6, 6);

            trafo.zero();
            inv_trafo.zero();

            trafo(0, 0) = 1;
            trafo(0, 3) = alpha;
            trafo(0, 5) = alpha;

            trafo(1, 1) = 1;
            trafo(1, 3) = alpha;
            trafo(1, 4) = alpha;

            trafo(2, 2) = 1;
            trafo(2, 4) = alpha;
            trafo(2, 5) = alpha;

            trafo(3, 3) = (1 - 2 * alpha);
            trafo(4, 4) = (1 - 2 * alpha);
            trafo(5, 5) = (1 - 2 * alpha);

            /////////////////////////////////////////////

            inv_trafo(0, 0) = 1;
            inv_trafo(0, 3) = (1 - 2 * alpha);
            inv_trafo(0, 5) = (1 - 2 * alpha);

            inv_trafo(1, 1) = 1;
            inv_trafo(1, 3) = (1 - 2 * alpha);
            inv_trafo(1, 4) = (1 - 2 * alpha);

            inv_trafo(2, 2) = 1;
            inv_trafo(2, 4) = (1 - 2 * alpha);
            inv_trafo(2, 5) = (1 - 2 * alpha);

            inv_trafo(3, 3) = 1. / (1 - 2 * alpha);
            inv_trafo(4, 4) = 1. / (1 - 2 * alpha);
            inv_trafo(5, 5) = 1. / (1 - 2 * alpha);
            return true;
        }

        if (el_type == libMesh::QUAD8) {
            trafo.resize(8, 8);
            inv_trafo.resize(8, 8);

            trafo.zero();
            inv_trafo.zero();

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

            trafo(4, 4) = (1 - 2 * alpha);
            trafo(5, 5) = (1 - 2 * alpha);
            trafo(6, 6) = (1 - 2 * alpha);
            trafo(7, 7) = (1 - 2 * alpha);

            /////////////////////////////////////////////

            inv_trafo(0, 0) = 1;
            inv_trafo(0, 4) = (1 - 2 * alpha);
            inv_trafo(0, 7) = (1 - 2 * alpha);

            inv_trafo(1, 1) = 1;
            inv_trafo(1, 4) = (1 - 2 * alpha);
            inv_trafo(1, 5) = (1 - 2 * alpha);

            inv_trafo(2, 2) = 1;
            inv_trafo(2, 5) = (1 - 2 * alpha);
            inv_trafo(2, 6) = (1 - 2 * alpha);

            inv_trafo(3, 3) = 1;
            inv_trafo(3, 6) = (1 - 2 * alpha);
            inv_trafo(3, 7) = (1 - 2 * alpha);

            inv_trafo(4, 4) = 1. / (1 - 2 * alpha);
            inv_trafo(5, 5) = 1. / (1 - 2 * alpha);
            inv_trafo(6, 6) = 1. / (1 - 2 * alpha);
            inv_trafo(7, 7) = 1. / (1 - 2 * alpha);
            return true;
        }

        if (el_type == libMesh::TET10) {
            trafo.resize(10, 10);
            inv_trafo.resize(10, 10);

            trafo.zero();
            inv_trafo.zero();

            trafo(0, 0) = 1;
            trafo(0, 4) = alpha;
            trafo(0, 6) = alpha;
            trafo(0, 7) = alpha;

            trafo(1, 1) = 1;
            trafo(1, 4) = alpha;
            trafo(1, 5) = alpha;
            trafo(1, 8) = alpha;

            trafo(2, 2) = 1;
            trafo(2, 5) = alpha;
            trafo(2, 6) = alpha;
            trafo(2, 9) = alpha;

            trafo(3, 3) = 1;
            trafo(3, 7) = alpha;
            trafo(3, 8) = alpha;
            trafo(3, 9) = alpha;

            trafo(4, 4) = (1 - 2 * alpha);
            trafo(5, 5) = (1 - 2 * alpha);
            trafo(6, 6) = (1 - 2 * alpha);
            trafo(7, 7) = (1 - 2 * alpha);
            trafo(8, 8) = (1 - 2 * alpha);
            trafo(9, 9) = (1 - 2 * alpha);

            /////////////////////////////////////////////

            inv_trafo(0, 0) = 1;
            inv_trafo(0, 4) = (1 - 2 * alpha);
            inv_trafo(0, 6) = (1 - 2 * alpha);
            inv_trafo(0, 7) = (1 - 2 * alpha);

            inv_trafo(1, 1) = 1;
            inv_trafo(1, 4) = (1 - 2 * alpha);
            inv_trafo(1, 5) = (1 - 2 * alpha);
            inv_trafo(1, 8) = (1 - 2 * alpha);

            inv_trafo(2, 2) = 1;
            inv_trafo(2, 5) = (1 - 2 * alpha);
            inv_trafo(2, 6) = (1 - 2 * alpha);
            inv_trafo(2, 9) = (1 - 2 * alpha);

            inv_trafo(3, 3) = 1;
            inv_trafo(3, 7) = (1 - 2 * alpha);
            inv_trafo(3, 8) = (1 - 2 * alpha);
            inv_trafo(3, 9) = (1 - 2 * alpha);

            inv_trafo(4, 4) = 1. / (1 - 2 * alpha);
            inv_trafo(5, 5) = 1. / (1 - 2 * alpha);
            inv_trafo(6, 6) = 1. / (1 - 2 * alpha);
            inv_trafo(7, 7) = 1. / (1 - 2 * alpha);
            inv_trafo(8, 8) = 1. / (1 - 2 * alpha);
            inv_trafo(9, 9) = 1. / (1 - 2 * alpha);
            return true;
        }

        assert(false && "not implemented");
        return false;
    }

    void DualBasis::assemble_biorth_weights(const libMesh::Elem &el,
                                            const int el_order,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize) {
        const auto dim = el.dim();
        std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(dim, libMesh::Order(el_order));

        const int order = order_for_l2_integral(dim, el, el_order, el, el_order);
        libMesh::QGauss qg(dim, libMesh::Order(order));
        biorth_elem->attach_quadrature_rule(&qg);
        biorth_elem->reinit(&el);
        assemble_biorth_weights(*biorth_elem, weights, normalize);
    }

    void DualBasis::assemble_biorth_weights(const libMesh::FEBase &fe,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize) {
        const auto &test = fe.get_phi();
        const auto &JxW = fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        libMesh::DenseMatrix<libMesh::Real> elmat;
        elmat.resize(n_test, n_test);
        elmat.zero();

        for (uint i = 0; i < n_test; ++i) {
            for (uint j = 0; j < n_test; ++j) {
                for (uint qp = 0; qp < n_qp; ++qp) {
                    elmat(i, j) += test[i][qp] * test[j][qp] * JxW[qp];
                }
            }
        }

        assemble_biorth_weights(elmat, weights, normalize);
    }

    void DualBasis::assemble_biorth_weights(const libMesh::Elem &el,
                                            const int el_order,
                                            const libMesh::DenseMatrix<libMesh::Real> &trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize) {
        const auto dim = el.dim();
        std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(dim, libMesh::Order(el_order));

        const int order = order_for_l2_integral(dim, el, el_order, el, el_order);
        libMesh::QGauss qg(dim, libMesh::Order(order));
        biorth_elem->attach_quadrature_rule(&qg);
        biorth_elem->reinit(&el);
        assemble_biorth_weights(*biorth_elem, trafo, weights, normalize);
    }

    void DualBasis::assemble_biorth_weights(const libMesh::FEBase &fe,
                                            const libMesh::DenseMatrix<libMesh::Real> &trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize) {
        const auto &test = fe.get_phi();
        const auto &JxW = fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        libMesh::DenseMatrix<libMesh::Real> elmat;
        elmat.resize(n_test, n_test);
        elmat.zero();

        std::vector<std::vector<double>> trafo_phi(n_test);
        ;
        for (uint i = 0; i < n_test; ++i) {
            trafo_phi[i].resize(n_qp);

            for (uint qp = 0; qp < n_qp; ++qp) {
                auto &val = trafo_phi[i][qp];
                val = 0.;

                for (uint j = 0; j < n_test; ++j) {
                    val += trafo(i, j) * test[j][qp];
                }
            }
        }

        for (uint i = 0; i < n_test; ++i) {
            for (uint j = 0; j < n_test; ++j) {
                for (uint qp = 0; qp < n_qp; ++qp) {
                    elmat(i, j) += trafo_phi[i][qp] * trafo_phi[j][qp] * JxW[qp];
                }
            }
        }

        assemble_biorth_weights(elmat, weights, normalize);
    }

    void DualBasis::assemble_biorth_weights(libMesh::DenseMatrix<libMesh::Real> &elmat,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize) {
        auto n_test = elmat.n();
        libMesh::DenseVector<libMesh::Real> rhs(n_test);
        rhs.zero();

        libMesh::DenseVector<libMesh::Real> sol(n_test);
        sol.zero();

        libMesh::DenseVector<libMesh::Real> sum_elmat(n_test);
        sum_elmat.zero();

        weights.resize(n_test, n_test);
        weights.zero();

        for (uint i = 0; i < n_test; ++i) {
            for (uint j = 0; j < n_test; ++j) {
                sum_elmat(i) += elmat(i, j);
            }
        }

        for (uint i = 0; i < n_test; ++i) {
            if (std::abs(sum_elmat(i)) < 1e-16) {
                sum_elmat(i) = 0;

                // set identity row where not defined
                for (uint j = 0; j < n_test; ++j) {
                    elmat(i, j) = (i == j);
                }
            }
        }

        for (uint i = 0; i < n_test; ++i) {
            if (sum_elmat(i) == 0) {
                continue;
            }

            rhs(i) = sum_elmat(i);

            elmat.cholesky_solve(rhs, sol);

            for (uint j = 0; j < n_test; ++j) {
                weights(i, j) = sol(j);
            }

            rhs(i) = 0;
        }

        if (normalize) {
            // normalization for consistently scaled coefficients
            for (uint i = 0; i < n_test; ++i) {
                if (sum_elmat(i) == 0) {
                    continue;
                }

                libMesh::Real t = 0;
                for (uint j = 0; j < n_test; ++j) {
                    t += weights(i, j);
                }

                for (uint j = 0; j < n_test; ++j) {
                    weights(i, j) *= 1. / t;
                }
            }
        }
    }

}  // namespace utopia
