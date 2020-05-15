#ifndef UTOPIA_LOCAL_L2_ASSEMBLER_R_HPP
#define UTOPIA_LOCAL_L2_ASSEMBLER_R_HPP

#include "MortarAssemble.hpp"
#include "utopia_LocalAssembler.hpp"
#include "utopia_QuadratureBasedAssembler.hpp"

#include <vector>

namespace utopia {
    class QMortarBuilder;

    class BidirectionalL2LocalAssembler final : public LocalAssembler, public QuadratureBasedAssembler {
    public:
        using Matrix = LocalAssembler::Matrix;

        BidirectionalL2LocalAssembler(const int dim, const bool use_biorth, const bool assemble_mass_mat = false);
        ~BidirectionalL2LocalAssembler();

        /**
         * @brief if you are performing volume to surface transfer
         * the method does not provide reliable results if the volume
         * element has facets aligned with the surface ones
         */
        bool assemble(const Elem &trial, FEType trial_type, const Elem &test, FEType test_type, Matrix &mats) override;

        bool assemble(const Elem &trial,
                      FEType trial_type,
                      const Elem &test,
                      FEType test_type,
                      std::vector<Matrix> &mat) override;

        inline const QMortarBuilder &get_q_builder() const {
            assert(q_builder);
            return *q_builder;
        }

        inline int n_forms() const override { return (assemble_mass_mat_) ? 4 : 2; }

        inline Type type(const int index) const override {
            assert(index < n_forms());
            assert(index >= 0);

            if (assemble_mass_mat_) {
                switch (index) {
                    case 0: {
                        return MASTER_X_SLAVE;
                    }
                    case 1: {
                        return SLAVE_X_SLAVE;
                    }
                    case 2: {
                        return SLAVE_X_MASTER;
                    }
                    case 3: {
                        return MASTER_X_MASTER;
                    }
                    default: {
                        assert(false);
                        return MASTER_X_SLAVE;
                    }
                }
            } else {
                switch (index) {
                    case 0: {
                        return MASTER_X_SLAVE;
                    }
                    case 1: {
                        return SLAVE_X_MASTER;
                    }
                    default: {
                        assert(false);
                        return MASTER_X_SLAVE;
                    }
                }
            }
        }

        void print_stats(std::ostream &os = std::cout) const override;

    private:
        int dim;
        bool use_biorth;
        bool must_compute_biorth;
        QMortar composite_ir;
        QMortar q_trial;
        QMortar q_test;

        Matrix test_biorth_weights, trial_biorth_weights;

        // std::shared_ptr<QMortarBuilder> q_builder;
        std::unique_ptr<libMesh::FEBase> trial_fe, test_fe;

        bool assemble_mass_mat_;

        void init_fe(const Elem &trial, FEType trial_type, const Elem &test, FEType test_type);

        void init_biorth(const Elem &trial, FEType trial_type, const Elem &test, FEType test_type);

        static void assemble_biorth_weights(const libMesh::Elem &el,
                                            const int dim,
                                            const libMesh::FEType &var_type,
                                            const int el_order,
                                            libMesh::DenseMatrix<libMesh::Real> &weights);
    };
}  // namespace utopia

#endif  // UTOPIA_LOCAL_L2_ASSEMBLER_HPP
