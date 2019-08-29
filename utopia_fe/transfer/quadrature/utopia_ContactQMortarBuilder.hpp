// #ifndef UTOPIA_CONTACT_Q_MORTAR_BUILDER_HPP
// #define UTOPIA_CONTACT_Q_MORTAR_BUILDER_HPP

// #include "libmesh/elem.h"
// #include "libmesh/fe.h"
// #include "libmesh/dense_matrix.h"
// #include "MortarAssemble.hpp"
// #include "utopia_Intersect.hpp"
// #include "utopia_Polygon2.hpp"

// #include <cassert>

// namespace utopia {

//     class ContactQMortarBuilder {
//     public:
//         using Elem = libMesh::Elem;
//         using FEType = libMesh::FEType;
//         using Real = libMesh::Real;
//         using Point = libMesh::Point;

//         virtual ~ContactQMortarBuilder() {}

//         virtual bool build(
//             const Elem &trial,
//             FEType trial_type,
//             const int trial_side_num,
//             const Elem &test,
//             FEType test_type,
//             const int test_side_num,
//             QMortar &q_trial,
//             QMortar &q_test) = 0;

//         virtual double get_total_intersection_volume() const = 0;
//         virtual bool is_affine() const = 0;
//         virtual const std::vector<double> &gap() const = 0;
//         virtual const std::vector<Point> &normal() const = 0;
//     };

//     class AffineContactQMortarBuilder3 final : public ContactQMortarBuilder {
//     public:
//         using Elem = libMesh::Elem;
//         using FEType = libMesh::FEType;
//         using Real = libMesh::Real;
//         using Point = libMesh::Point;
//         using Matrix = libMesh::DenseMatrix<libMesh::Real>;

//         ~AffineContactQMortarBuilder3() {}
//         AffineContactQMortarBuilder3(const Real search_radius);

//         bool build(
//             const Elem &trial,
//             FEType trial_type,
//             const int trial_side_num,
//             const Elem &test,
//             FEType test_type,
//             const int test_side_num,
//             QMortar &q_trial,
//             QMortar &q_test) override;

//         double get_total_intersection_volume() const override;
//         inline bool is_affine() const override { return true; }

//         inline const std::vector<double> &gap() const override
//         {
//             return gap_;
//         }

//         inline const std::vector<Point> &normal() const override
//         {
//             return normal_;
//         }

//     private:
//         Matrix trial_polygon, test_polygon;
//         Matrix trial_isect, test_isect;
//         QMortar trial_ir, test_ir;

//         Point trial_n, test_n;

//         Box trial_box, test_box;

//         Real total_intersection_volume;
//         Real search_radius;

//         std::vector<double> gap_;
//         std::vector<Point> normal_;
//     };

//     class WarpedContactQMortarBuilder3 final : public ContactQMortarBuilder {
//     public:
//         using Elem = libMesh::Elem;
//         using FEType = libMesh::FEType;
//         using Real = libMesh::Real;
//         using Point = libMesh::Point;
//         using Matrix = libMesh::DenseMatrix<libMesh::Real>;

//         ~WarpedContactQMortarBuilder3() {}
//         WarpedContactQMortarBuilder3(const Real search_radius);

//         bool build(
//             const Elem &trial,
//             FEType trial_type,
//             const int trial_side_num,
//             const Elem &test,
//             FEType test_type,
//             const int test_side_num,
//             QMortar &q_trial,
//             QMortar &q_test) override;

//         double get_total_intersection_volume() const override;
//         inline bool is_affine() const override { return false; }

//         inline const std::vector<double> &gap() const override
//         {
//             return gap_;
//         }

//         inline const std::vector<Point> &normal() const override
//         {
//             return normal_;
//         }

//     private:
//         Matrix trial_polygon, test_polygon;
//         Matrix trial_isect, test_isect;
//         QMortar trial_ir, test_ir;

//         Point trial_n, test_n;

//         Box trial_box, test_box;

//         Real total_intersection_volume;
//         Real search_radius;

//         Polygon3 trial_poly3_, test_poly3_;

//         int ref_quad_order_;
//         std::vector<Polygon2::Vector> tri_q_points_;
//         std::vector<Polygon2::Scalar> tri_q_weights_;

//         std::vector<Polygon3::Vector> composite_q_points_;
//         std::vector<Polygon3::Scalar> composite_q_weights_;
//         std::vector<double> trial_gap_, test_gap_, gap_;
//         std::vector<Point> normal_;

//         void init_ref_quad(const int order);
//     };

// }

// #endif //UTOPIA_CONTACT_Q_MORTAR_BUILDER_HPP