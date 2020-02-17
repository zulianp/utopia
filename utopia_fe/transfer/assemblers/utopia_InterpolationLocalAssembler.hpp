#ifndef UTOPIA_INTERPOLATION_LOCAL_ASSEMBLER_HPP
#define UTOPIA_INTERPOLATION_LOCAL_ASSEMBLER_HPP

#include "utopia_LocalAssembler.hpp"
#include "MortarAssemble.hpp"

#include "libmesh/dense_matrix.h"
#include "libmesh/point.h"


namespace utopia {
    class InterpolationLocalAssembler final : public LocalAssembler {
    public:
        using Point = libMesh::Point;
        using Matrix = libMesh::DenseMatrix<libMesh::Real>;

        InterpolationLocalAssembler(const int dim, const bool nested_meshes = false)
        : dim(dim), nested_meshes(nested_meshes), tol(1e-10), force_affine_(false)
        { }

        bool assemble(
            const Elem &trial,
            FEType trial_type,
            const Elem &test,
            FEType test_type,
            Matrix &mat
            ) override;

        void force_affine(const bool val) {
            force_affine_ = val;
        }
        
    private:
        int dim;
        bool nested_meshes;
        double tol;
        bool force_affine_;
        std::shared_ptr<QMortar> q_trial, q_test;
        Matrix trial_pts, test_pts;

        std::shared_ptr<Transform> get_trafo(const Elem &elem) const;
        bool check_valid(const Matrix &mat) const;
        void init_q(const std::size_t n_qp);
        void contained_points(const Elem &trial, const Elem &test, std::vector<int> &test_dofs);
        void contained_points_1(const Elem &trial, const Elem &test, std::vector<int> &test_dofs);
        void contained_points_2(const Elem &trial, const Elem &test, std::vector<int> &test_dofs);
        void contained_points_3(const Elem &trial, const Elem &test, std::vector<int> &test_dofs) const;

    };
}

#endif //UTOPIA_INTERPOLATION_LOCAL_ASSEMBLER_HPP
