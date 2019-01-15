#ifndef UTOPIA_Q_MORTAR_BUILDER_HPP
#define UTOPIA_Q_MORTAR_BUILDER_HPP

#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/dense_matrix.h"
#include "MortarAssemble.hpp"

#include <cassert>

namespace utopia {

	class QMortarBuilder {
	public:
		using Elem = libMesh::Elem;
		using FEType = libMesh::FEType;

		virtual ~QMortarBuilder() {}

		virtual bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) = 0;

		virtual double get_total_intersection_volume() const = 0;
	};


	class QMortarBuilder1 final : public QMortarBuilder {
	public:
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;
		using Real = libMesh::Real;

		QMortarBuilder1()
		: total_intersection_volume(0), composite_ir(1)
		{}

		bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) override;

		inline double get_total_intersection_volume() const override
		{
			return total_intersection_volume;
		}

	private:
		libMesh::Point intersection[2];
		libMesh::Point u, v, w, min_p, max_p, min_q, max_q, r;

		Real total_intersection_volume;
		QMortar composite_ir;
	};


	class QMortarBuilder2 final : public QMortarBuilder {
	public:
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;
		using Real = libMesh::Real;

		QMortarBuilder2()
		: total_intersection_volume(0), composite_ir(2)
		{}

		bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) override;

		inline double get_total_intersection_volume() const override
		{
			return total_intersection_volume;
		}

	private:
		Matrix trial_pts;
		Matrix test_pts;
		Matrix intersection;

		Real total_intersection_volume;
		QMortar composite_ir;
	};

	class QMortarBuilderShell2 final : public QMortarBuilder {
	public:
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;
		using Real = libMesh::Real;

		QMortarBuilderShell2()
		: total_intersection_volume(0), composite_ir(2)
		{}

		bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) override;

		inline double get_total_intersection_volume() const override
		{
			return total_intersection_volume;
		}

	private:
		Matrix trial_pts;
		Matrix test_pts;
		Matrix intersection;

		Real total_intersection_volume;
		QMortar composite_ir;

		libMesh::Point trial_normal, test_normal;

		//intersector
		typedef Intersector::Scalar Scalar;
		
		Scalar A   [3 * 3], b   [3];
		Scalar Ainv[3 * 3], binv[3];

		Matrix ref_trial_pts;
		Matrix ref_test_pts;
		Matrix ref_trial_pts_2;
		Matrix ref_test_pts_2;
		Matrix ref_intersection_2;

		Matrix ref_intersection_slave;
		Matrix ref_intersection_master;

		inline void compute_normal(const Elem &elem, libMesh::Point &n) const
		{
			using namespace libMesh;
			Point o, u, v;

			o = elem.point(0);
			u = elem.point(1);
			v = elem.point(2);
			u -= o;
			v -= o;
			n = u.cross(v);

			n *= 1./n.norm();
		}

		bool intersect();
	};

	class QMortarBuilder3 final : public QMortarBuilder {
	public:
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;
		using Real = libMesh::Real;

		QMortarBuilder3()
		: total_intersection_volume(0.), composite_ir(3)
		{}

		bool build(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test) override;

		inline double get_total_intersection_volume() const override
		{
			return total_intersection_volume;
		}

	private:
		Real total_intersection_volume;
		QMortar composite_ir;

		Polyhedron trial_poly, test_poly;
		Polyhedron intersection, temp_poly;

		Matrix shell_poly;

		bool build_vol_2_surf(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			QMortar &q_trial,
			QMortar &q_test);
	};
}


#endif //UTOPIA_Q_MORTAR_BUILDER_HPP
