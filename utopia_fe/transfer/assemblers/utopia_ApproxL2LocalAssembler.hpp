#ifndef UTOPIA_APPROX_L2_LOCAL_ASSEMBLER_HPP
#define UTOPIA_APPROX_L2_LOCAL_ASSEMBLER_HPP

#include "utopia_LocalAssembler.hpp"
#include "MortarAssemble.hpp"

#include "libmesh/dense_matrix.h"
#include "libmesh/point.h"


namespace utopia {
    class ApproxL2LocalAssembler final : public LocalAssembler {
    public:
        using Point = libMesh::Point;
        using Matrix = libMesh::DenseMatrix<libMesh::Real>;

        ApproxL2LocalAssembler(
            const int dim,
            const bool nested_meshes = false)
        : dim(dim), quadrature_order(-1), nested_meshes(nested_meshes), tol(1e-10)
        { }

        bool assemble(
            const Elem &trial,
            FEType trial_type,
            const Elem &test,
            FEType test_type,
            Matrix &mat
            ) override;

        bool assemble(
                const Elem &trial,
                FEType trial_type,
                const Elem &test,
                FEType test_type,
                std::vector<Matrix> &mat
                ) override;

        inline void set_quadrature_order(const int val)
        {
            quadrature_order = val;
        }

        inline int n_forms() const override
        {
            return 2;
        }

        inline Type type(const int index) const override
        {
            assert(index < n_forms());
            assert(index >= 0);

            return index == 0 ? MASTER_X_SLAVE : SLAVE_X_SLAVE;
        }

    private:
        int dim;
        int quadrature_order;
        bool nested_meshes;
        double tol;
        std::shared_ptr<QMortar> q_trial, q_test;
        Matrix trial_pts;

        std::shared_ptr<Transform> get_trafo(const Elem &elem) const;
        bool check_valid(const Matrix &mat) const;
        bool init_q(
        const libMesh::Elem &trial,
        libMesh::FEType trial_type,
        const libMesh::Elem &test,
        libMesh::FEType test_type);

        void contained_points(const libMesh::Elem &trial,   const libMesh::QBase &q, std::vector<int> &index);
        void contained_points_2(const libMesh::Elem &trial, const libMesh::QBase &q, std::vector<int> &index);
        void contained_points_3(const libMesh::Elem &trial, const libMesh::QBase &q, std::vector<int> &index) const;

        void assemble(libMesh::FEBase &trial, libMesh::FEBase &test, Matrix &mat) const;
    };
}


#endif //UTOPIA_APPROX_L2_LOCAL_ASSEMBLER_HPP
