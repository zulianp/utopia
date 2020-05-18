#include "utopia_libmesh_Transform.hpp"
#include <libmesh/fe.h>
#include "utopia_intersector.hpp"

namespace utopia {

    void Transform1::apply(const libMesh::Point &ref, libMesh::Point &world) const {
        assert((libMesh::FE<1, libMesh::LAGRANGE>::on_reference_element(ref, elem_.type(), 1e-3)));
        world = libMesh::FE<1, libMesh::LAGRANGE>::map(&elem_, ref);
    }

    void Transform1::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const {
        ref = libMesh::FE<1, libMesh::LAGRANGE>::inverse_map(&elem_, world, 1e-10);
        assert((libMesh::FE<1, libMesh::LAGRANGE>::on_reference_element(ref, elem_.type(), 1e-3)));
        assert((libMesh::FE<1, libMesh::LAGRANGE>::map(&elem_, ref).absolute_fuzzy_equals(world, 1e-8)));
    }

    void Transform2::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const {
        ref = libMesh::FE<2, libMesh::LAGRANGE>::inverse_map(&elem_, world, 1e-10);
        assert((libMesh::FE<2, libMesh::LAGRANGE>::on_reference_element(ref, elem_.type(), 1e-3)));
        assert((libMesh::FE<2, libMesh::LAGRANGE>::map(&elem_, ref).absolute_fuzzy_equals(world, 1e-8)));
    }

    void Transform2::apply(const libMesh::Point &ref, libMesh::Point &world) const {
        world = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem_, ref);
    }

    void AffineTransform2::compute_affine_transformation(const libMesh::Elem &elem,
                                                         libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                         libMesh::DenseVector<libMesh::Real> &A_inv_m_b) {
        libMesh::Point p0, p1, p2;

        A_inv.resize(2, 2);
        A_inv_m_b.resize(2);

        libMesh::DenseMatrix<libMesh::Real> A;

        A.resize(2, 2);

        std::vector<const libMesh::Node *> elem_nodes;

        std::vector<libMesh::Point> reference_points;

        libMesh::Point ref_p0(0.0, 0.0, 0.0);

        libMesh::Point ref_p1(1.0, 0.0, 0.0);

        libMesh::Point ref_p2(0.0, 1.0, 0.0);

        p0 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, ref_p0);

        p1 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, ref_p1);

        p2 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, ref_p2);

        A(0, 0) = (p1(0) - p0(0));
        A(0, 1) = (p2(0) - p0(0));
        A(1, 0) = (p1(1) - p0(1));
        A(1, 1) = (p2(1) - p0(1));

        libMesh::Real det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

        A_inv(0, 0) = 1. / det * A(1, 1);
        A_inv(1, 1) = 1. / det * A(0, 0);
        A_inv(0, 1) = -1. / det * A(0, 1);
        A_inv(1, 0) = -1. / det * A(1, 0);

        A_inv_m_b(0) = -1.0 * A_inv(0, 0) * p0(0) - A_inv(0, 1) * p0(1);
        A_inv_m_b(1) = -1.0 * A_inv(1, 0) * p0(0) - A_inv(1, 1) * p0(1);
    }

    void AffineTransform2::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const {
        ref(0) = A_inv_m_b_(0) + A_inv_(0, 0) * world(0) + A_inv_(0, 1) * world(1);
        ref(1) = A_inv_m_b_(1) + A_inv_(1, 0) * world(0) + A_inv_(1, 1) * world(1);
        ref(2) = 0.0;
    }

    void Transform3::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const {
        ref = libMesh::FE<3, libMesh::LAGRANGE>::inverse_map(&elem_, world);
        // assert( (libMesh::FE<3, libMesh::LAGRANGE>::on_reference_element(ref, elem_.type(), 1e-6)) );
        assert((libMesh::FE<3, libMesh::LAGRANGE>::map(&elem_, ref).absolute_fuzzy_equals(world, 1e-8)));
    }

    void Transform3::apply(const libMesh::Point &ref, libMesh::Point &world) const {
        world = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem_, ref);
    }

    void AffineTransform3::compute_affine_transformation(const libMesh::Elem &elem,
                                                         libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                         libMesh::DenseVector<libMesh::Real> &A_inv_m_b) {
        libMesh::Point p0, p1, p2, p3;

        libMesh::Point ref_p0(0.0, 0.0, 0.0);
        libMesh::Point ref_p1(1.0, 0.0, 0.0);
        libMesh::Point ref_p2(0.0, 1.0, 0.0);
        libMesh::Point ref_p3(0.0, 0.0, 1.0);

        A_inv.resize(3, 3);
        A_inv_m_b.resize(3);

        libMesh::DenseMatrix<libMesh::Real> A;

        A.resize(3, 3);

        std::vector<const libMesh::Node *> elem_nodes;

        std::vector<libMesh::Point> reference_points;

        p0 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p0);
        p1 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p1);
        p2 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p2);
        p3 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p3);

        A(0, 0) = p1(0) - p0(0);
        A(0, 1) = p2(0) - p0(0);
        A(0, 2) = p3(0) - p0(0);
        A(1, 0) = p1(1) - p0(1);
        A(1, 1) = p2(1) - p0(1);
        A(1, 2) = p3(1) - p0(1);
        A(2, 0) = p1(2) - p0(2);
        A(2, 1) = p2(2) - p0(2);
        A(2, 2) = p3(2) - p0(2);

        libMesh::Real det = Intersector::det_3(&A.get_values()[0]);
        Intersector::inverse_3(&A.get_values()[0], det, &A_inv.get_values()[0]);

        A_inv_m_b(0) = -1.0 * A_inv(0, 0) * p0(0) - A_inv(0, 1) * p0(1) - 1.0 * A_inv(0, 2) * p0(2);
        A_inv_m_b(1) = -1.0 * A_inv(1, 0) * p0(0) - A_inv(1, 1) * p0(1) - 1.0 * A_inv(1, 2) * p0(2);
        A_inv_m_b(2) = -1.0 * A_inv(2, 0) * p0(0) - A_inv(2, 1) * p0(1) - 1.0 * A_inv(2, 2) * p0(2);
    }

    void SideAffineTransform2::compute_affine_transformation(const libMesh::Elem &elem,
                                                             const int side,
                                                             libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                             libMesh::DenseVector<libMesh::Real> &A_inv_m_b) {
        // ref element -1, 1
        libMesh::Point u, n;
        auto side_ptr = elem.build_side_ptr(side);

        libMesh::Point ref_p0(-1.);
        libMesh::Point ref_p1(1.0);

        libMesh::Point p0 = libMesh::FE<1, libMesh::LAGRANGE>::map(side_ptr.get(), ref_p0);
        libMesh::Point p1 = libMesh::FE<1, libMesh::LAGRANGE>::map(side_ptr.get(), ref_p1);

        u = p1 - p0;

        n(0) = -u(1);
        n(1) = u(0);

        n /= n.norm();

        A_inv.resize(2, 2);
        A_inv_m_b.resize(2);

        libMesh::DenseMatrix<libMesh::Real> A;
        A.resize(2, 2);

        A(0, 0) = u(0);
        A(0, 1) = n(0);
        A(1, 0) = u(1);
        A(1, 1) = n(1);

        const libMesh::Real det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

        A_inv(0, 0) = 1. / det * A(1, 1);
        A_inv(1, 1) = 1. / det * A(0, 0);
        A_inv(0, 1) = -1. / det * A(0, 1);
        A_inv(1, 0) = -1. / det * A(1, 0);

        A_inv_m_b(0) = -1.0 * A_inv(0, 0) * p0(0) - A_inv(0, 1) * p0(1);
        A_inv_m_b(1) = -1.0 * A_inv(1, 0) * p0(0) - A_inv(1, 1) * p0(1);
    }

    void SideAffineTransform3::compute_affine_transformation(const libMesh::Elem &elem,
                                                             const int side,
                                                             libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                             libMesh::DenseVector<libMesh::Real> &A_inv_m_b) {
        auto side_ptr = elem.build_side_ptr(side);

        libMesh::Point ref_p0(0.0, 0.0);
        libMesh::Point ref_p1(1.0, 0.0);
        libMesh::Point ref_p2(0.0, 1.0);

        libMesh::Point p0 = libMesh::FE<2, libMesh::LAGRANGE>::map(side_ptr.get(), ref_p0);
        libMesh::Point p1 = libMesh::FE<2, libMesh::LAGRANGE>::map(side_ptr.get(), ref_p1);
        libMesh::Point p2 = libMesh::FE<2, libMesh::LAGRANGE>::map(side_ptr.get(), ref_p2);

        libMesh::Point u, v, n;

        u = p1 - p0;
        v = p2 - p0;

        n = u.cross(v);
        n /= n.norm();

        A_inv.resize(3, 3);
        A_inv_m_b.resize(3);

        libMesh::DenseMatrix<libMesh::Real> A;

        A.resize(3, 3);

        A(0, 0) = u(0);
        A(0, 1) = v(0);
        A(0, 2) = n(0);
        A(1, 0) = u(1);
        A(1, 1) = v(1);
        A(1, 2) = n(1);
        A(2, 0) = u(2);
        A(2, 1) = v(2);
        A(2, 2) = n(2);

        libMesh::Real det = Intersector::det_3(&A.get_values()[0]);
        Intersector::inverse_3(&A.get_values()[0], det, &A_inv.get_values()[0]);

        A_inv_m_b(0) = -1.0 * A_inv(0, 0) * p0(0) - A_inv(0, 1) * p0(1) - 1.0 * A_inv(0, 2) * p0(2);
        A_inv_m_b(1) = -1.0 * A_inv(1, 0) * p0(0) - A_inv(1, 1) * p0(1) - 1.0 * A_inv(1, 2) * p0(2);
        A_inv_m_b(2) = -1.0 * A_inv(2, 0) * p0(0) - A_inv(2, 1) * p0(1) - 1.0 * A_inv(2, 2) * p0(2);
    }

    void AffineTransform3::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const {
        ref(0) = A_inv_m_b_(0) + A_inv_(0, 0) * world(0) + A_inv_(0, 1) * world(1) + A_inv_(0, 2) * world(2);
        ref(1) = A_inv_m_b_(1) + A_inv_(1, 0) * world(0) + A_inv_(1, 1) * world(1) + A_inv_(1, 2) * world(2);
        ref(2) = A_inv_m_b_(2) + A_inv_(2, 0) * world(0) + A_inv_(2, 1) * world(1) + A_inv_(2, 2) * world(2);
    }

}  // namespace utopia
