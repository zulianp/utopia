#ifndef UTOPIA_DUAL_BASIS_HPP
#define UTOPIA_DUAL_BASIS_HPP

#include "utopia_fe_base.hpp"

#include "libmesh/dense_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/reference_elem.h"

namespace utopia {

    bool is_diag(const libMesh::DenseMatrix<double> &d, const bool verbose = true);

    //@brief from the paper DUAL QUADRATIC MORTAR FINITE ELEMENT METHODS FOR 3D FINITE DEFORMATION CONTACT∗
    class DualBasis {
    public:
        using VectorValueT = libMesh::VectorValue<double>;

        libMesh::DenseMatrix<libMesh::Real> trafo_;
        libMesh::DenseMatrix<libMesh::Real> inv_trafo_;
        libMesh::DenseMatrix<libMesh::Real> weights_;

        std::vector<std::vector<double>> phi_;
        std::vector<std::vector<VectorValueT>> dphi_;

        int order = -1;

        bool compute_phi = false;
        bool compute_dphi = false;
        bool compute_divphi = false;

        inline bool empty() const { return weights_.n() == 0; }

        inline bool must_compute_values() const { return compute_dphi || compute_dphi || compute_divphi; }

        constexpr static const double DEFAULT_ALPHA = 1. / 5;

        void init(const libMesh::ElemType &type, const double alpha = DEFAULT_ALPHA);

        void compute_values(const libMesh::FEBase &fe);

        static bool build_trafo_and_weights(const libMesh::ElemType type,
                                            const int order,
                                            const double alpha,
                                            libMesh::DenseMatrix<libMesh::Real> &trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &inv_trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &weights);

        template <class EDofMap>
        static bool build_global_trafo(const libMesh::MeshBase &mesh,
                                       const SizeType n_local_dofs,
                                       const EDofMap &dof_map,
                                       const UVector &elem_to_transform,
                                       const double alpha,
                                       USparseMatrix &mat,
                                       const bool inverse = false);
        // transposed trafo
        static bool assemble_local_trafo(const libMesh::ElemType el_type,
                                         const double alpha,
                                         libMesh::DenseMatrix<libMesh::Real> &trafo,
                                         libMesh::DenseMatrix<libMesh::Real> &inv_trafo);

        static void assemble_biorth_weights(const libMesh::Elem &el,
                                            const int el_order,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize = true);

        static void assemble_biorth_weights(const libMesh::FEBase &fe,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize = true);

        static void assemble_biorth_weights(const libMesh::Elem &el,
                                            const int el_order,
                                            const libMesh::DenseMatrix<libMesh::Real> &trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize = true);

        static void assemble_biorth_weights(const libMesh::FEBase &fe,
                                            const libMesh::DenseMatrix<libMesh::Real> &trafo,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize = true);

        static void assemble_biorth_weights(libMesh::DenseMatrix<libMesh::Real> &elmat,
                                            libMesh::DenseMatrix<libMesh::Real> &weights,
                                            const bool normalize);
    };
}  // namespace utopia

#endif  // UTOPIA_DUAL_BASIS_HPP
