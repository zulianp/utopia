#ifndef UTOPIA_Q_MORTAR_BUILDER_HPP
#define UTOPIA_Q_MORTAR_BUILDER_HPP

#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/dense_matrix.h"
#include "MortarAssemble.hpp"

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

	class QMortarBuilder3 final : public QMortarBuilder {
	public:
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;
		using Real = libMesh::Real;

		QMortarBuilder3()
		: total_intersection_volume(0), composite_ir(3)
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
	};
}


#endif //UTOPIA_Q_MORTAR_BUILDER_HPP
