#include "utopia_BidirectionalL2LocalAssembler.hpp"
#include "utopia_QMortarBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <sstream>

namespace utopia {

    BidirectionalL2LocalAssembler::BidirectionalL2LocalAssembler(const int dim,
                                                                 const bool use_biorth,
                                                                 const bool assemble_mass_mat)
        : dim(dim),
          use_biorth(use_biorth),
          must_compute_biorth(use_biorth),
          composite_ir(dim),
          q_trial(dim),
          q_test(dim),
          assemble_mass_mat_(assemble_mass_mat) {
        if (dim == 2) {
            q_builder = std::make_shared<QMortarBuilder2>();
        } else {
            assert(dim == 3);
            q_builder = std::make_shared<QMortarBuilder3>();
        }
    }

    BidirectionalL2LocalAssembler::~BidirectionalL2LocalAssembler() {}

    bool BidirectionalL2LocalAssembler::assemble(const Elem &trial,
                                                 FEType trial_type,
                                                 const Elem &test,
                                                 FEType test_type,
                                                 Matrix &mat) {
        assert(false);
        return false;
    }

    bool BidirectionalL2LocalAssembler::assemble(const Elem &trial,
                                                 FEType trial_type,
                                                 const Elem &test,
                                                 FEType test_type,
                                                 std::vector<Matrix> &mat) {
        mat.resize(n_forms());

        auto trial_fe = libMesh::FEBase::build(trial.dim(), trial_type);
        auto test_fe = libMesh::FEBase::build(test.dim(), test_type);

        const int order = std::max(std::max(order_for_l2_integral(dim, trial, trial_type.order, test, test_type.order),
                                            order_for_l2_integral(dim, test, test_type.order, test, test_type.order)),
                                   order_for_l2_integral(dim, trial, trial_type.order, trial, trial_type.order));

        if (!q_builder->build(trial, trial_type, test, test_type, q_trial, q_test)) {
            return false;
        }

        init_biorth(trial, trial_type, test, test_type);

        init_fe(trial, trial_type, test, test_type);

        trial_fe->attach_quadrature_rule(&q_trial);
        trial_fe->get_phi();
        trial_fe->get_JxW();
        trial_fe->reinit(&trial);

        test_fe->attach_quadrature_rule(&q_test);
        test_fe->get_phi();
        test_fe->get_JxW();
        test_fe->reinit(&test);

        if (assemble_mass_mat_) {
            if (use_biorth) {
                mortar_assemble_weighted_biorth(*trial_fe, *test_fe, test_biorth_weights, mat[0]);
                mortar_assemble_weighted_biorth(*test_fe, *test_fe, test_biorth_weights, mat[1]);

                mortar_assemble_weighted_biorth(*test_fe, *trial_fe, trial_biorth_weights, mat[2]);
                mortar_assemble_weighted_biorth(*trial_fe, *trial_fe, trial_biorth_weights, mat[3]);
            } else {
                mortar_assemble(*trial_fe, *test_fe, mat[0]);
                mortar_assemble(*test_fe, *test_fe, mat[1]);

                mortar_assemble(*test_fe, *trial_fe, mat[2]);
                mortar_assemble(*trial_fe, *trial_fe, mat[3]);
            }

        } else {
            if (use_biorth) {
                mortar_assemble_weighted_biorth(*trial_fe, *test_fe, test_biorth_weights, mat[0]);
                mortar_assemble_weighted_biorth(*test_fe, *trial_fe, trial_biorth_weights, mat[1]);
            } else {
                mortar_assemble(*trial_fe, *test_fe, mat[0]);
                mortar_assemble(*test_fe, *trial_fe, mat[1]);
            }
        }

        return true;
    }

    void BidirectionalL2LocalAssembler::init_fe(const Elem &trial,
                                                FEType trial_type,
                                                const Elem &test,
                                                FEType test_type) {
        if (trial_fe) return;

        trial_fe = libMesh::FEBase::build(trial.dim(), trial_type);
        test_fe = libMesh::FEBase::build(test.dim(), test_type);
    }

    void BidirectionalL2LocalAssembler::init_biorth(const Elem &trial,
                                                    FEType trial_type,
                                                    const Elem &test,
                                                    FEType test_type) {
        if (!use_biorth) return;
        if (!must_compute_biorth) return;

        assemble_biorth_weights(test, test.dim(), test_type, test_type.order, test_biorth_weights);

        assemble_biorth_weights(trial, trial.dim(), trial_type, trial_type.order, trial_biorth_weights);

        must_compute_biorth = false;
    }

    void BidirectionalL2LocalAssembler::assemble_biorth_weights(const libMesh::Elem &el,
                                                                const int dim,
                                                                const libMesh::FEType &var_type,
                                                                const int el_order,
                                                                libMesh::DenseMatrix<libMesh::Real> &weights) {
        std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(dim, var_type);

        const int order = order_for_l2_integral(dim, el, el_order, el, el_order);

        libMesh::QGauss qg(dim, libMesh::Order(order));
        biorth_elem->attach_quadrature_rule(&qg);
        biorth_elem->reinit(&el);
        mortar_assemble_weights(*biorth_elem, weights);
    }

    void BidirectionalL2LocalAssembler::print_stats(std::ostream &os) const {}
}  // namespace utopia
