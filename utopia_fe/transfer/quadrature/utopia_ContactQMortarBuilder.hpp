#ifndef UTOPIA_CONTACT_Q_MORTAR_BUILDER_HPP
#define UTOPIA_CONTACT_Q_MORTAR_BUILDER_HPP

#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/dense_matrix.h"
#include "MortarAssemble.hpp"

#include <cassert>

namespace utopia {

	class ContactQMortarBuilder {
	public:
		using Elem = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Real = libMesh::Real;

		virtual ~ContactQMortarBuilder() {}

		virtual bool build(
			const Elem &trial,
			FEType trial_type,
			const int trial_side_num,
			const Elem &test,
			FEType test_type,
			const int test_side_num,
			QMortar &q_trial,
			QMortar &q_test) = 0;

		virtual double get_total_intersection_volume() const = 0;
	};

	class AffineContactQMortarBuilder3 final : public ContactQMortarBuilder {
	public:
		using Elem = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Real = libMesh::Real;
		using Point = libMesh::Point;
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;

		~AffineContactQMortarBuilder3() {}
		AffineContactQMortarBuilder3(const Real search_radius);

		bool build(
			const Elem &trial,
			FEType trial_type,
			const int trial_side_num,
			const Elem &test,
			FEType test_type,
			const int test_side_num,
			QMortar &q_trial,
			QMortar &q_test) override;

		double get_total_intersection_volume() const override;

	private:
		Matrix trial_polygon, test_polygon;
		Matrix trial_isect, test_isect;
		QMortar trial_ir, test_ir;

		Point trial_n, test_n;

		Box trial_box, test_box;

		Real total_intersection_volume; 
		Real search_radius;
	};

	class WarpedContactQMortarBuilder3 final : public ContactQMortarBuilder {
	public:
		using Elem = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Real = libMesh::Real;
		using Point = libMesh::Point;
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;

		~WarpedContactQMortarBuilder3() {}
		WarpedContactQMortarBuilder3(const Real search_radius);

		bool build(
			const Elem &trial,
			FEType trial_type,
			const int trial_side_num,
			const Elem &test,
			FEType test_type,
			const int test_side_num,
			QMortar &q_trial,
			QMortar &q_test) override;

		double get_total_intersection_volume() const override;

	private:
		Matrix trial_polygon, test_polygon;
		Matrix trial_isect, test_isect;
		QMortar trial_ir, test_ir;

		Point trial_n, test_n;

		Box trial_box, test_box;

		Real total_intersection_volume; 
		Real search_radius;
	};

}

#endif //UTOPIA_CONTACT_Q_MORTAR_BUILDER_HPP