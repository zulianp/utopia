#ifndef UTOPIA_MORTAR_ASSEMBLE_HPP
#define UTOPIA_MORTAR_ASSEMBLE_HPP



//#include "utopia_fe_core.hpp"
//#include "utopia_LibMeshBackend.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_LibMeshBackend.hpp"

#include "HashGrid.hpp"

#include <libmesh/elem.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/sparse_matrix.h>
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_quadrature_type.h"

// // Define the Finite Element object.
#include "libmesh/fe.h"

// // Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

#include "utopia_intersector.hpp"
#include "utopia_libmesh_Utils.hpp"

#include <memory>
#include <math.h>
#include <algorithm>

namespace utopia {

	void compute_side_normal(const int dim, const libMesh::Elem &side, libMesh::Point &n);

	int order_for_l2_integral(const int dim,
							  const libMesh::Elem &master_el,
							  const int master_order,
							  const libMesh::Elem &slave_el,
							  const int slave_order);

	class Transform {
	public:
		virtual ~Transform() {}
		virtual void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const = 0;
		virtual void apply(const libMesh::Point &ref, libMesh::Point &world) const = 0;
	};


	class Transform1 : public Transform {
	public:
		Transform1(const libMesh::Elem &elem)
		: elem_(elem)
		{}

		void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
		void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

	private:
		const libMesh::Elem &elem_;
	};


	class Transform2 : public Transform {
	public:
		// Transform2(const libMesh::DenseMatrix<libMesh::Real> &polygon)
		// : polygon_(polygon)
		// {
		// 	assert(polygon_.m() == 3 || polygon_.m() == 4 && "must be either a triangle or a quad");
		// }

		Transform2(const libMesh::Elem &elem)
		: elem_(elem)
		{}

		void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
		void apply(const libMesh::Point &ref, libMesh::Point &world) const override;

		static void compute_affine_transformation(const libMesh::Elem * elem, libMesh::DenseMatrix<libMesh::Real> &A_inv);

	private:

		const libMesh::Elem &elem_;

	};

	class Transform3 : public Transform {
	public:
	public:
		// Transform3(const Polyhedron &polyhedron, const libMesh::Elem &elem)
		// : polyhedron_(polyhedron), elem_(elem)
		// {
		// 	assert(polyhedron.n_nodes == 4 || polyhedron.n_nodes == 8 && "must be either a tetrahedron or a hex");
		// }

		Transform3(const libMesh::Elem &elem)
		: elem_(elem)
		{ }

		void transform_to_reference(const libMesh::Point &world, libMesh::Point &refm) const override;
		void apply(const libMesh::Point &ref, libMesh::Point &world) const override;


	private:
		// const Polyhedron &polyhedron_;
		const libMesh::Elem &elem_;

	};


	class AffineTransform2 : public Transform {
	public:
		AffineTransform2(const libMesh::Elem &elem)
		{
			compute_affine_transformation(elem, A_inv_ , A_inv_m_b_);
		}

		AffineTransform2(const libMesh::DenseMatrix<libMesh::Real> &A_inv,
						 const libMesh::DenseVector<libMesh::Real> A_inv_m_b)
		: A_inv_(A_inv), A_inv_m_b_(A_inv_m_b)
		{}

		void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
		void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

		AffineTransform2() {}


		libMesh::DenseMatrix<libMesh::Real> &A_inv()
		{
			return A_inv_;
		}

		libMesh::DenseVector<libMesh::Real> &A_inv_m_b()
		{
			return A_inv_m_b_;
		}

	private:
		libMesh::DenseMatrix<libMesh::Real> A_inv_;
		libMesh::DenseVector<libMesh::Real> A_inv_m_b_;



		static void compute_affine_transformation(const libMesh::Elem &elem,
												  libMesh::DenseMatrix<libMesh::Real> &A_inv,
												  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
	};


	class AffineTransform3 : public Transform {
	public:
		AffineTransform3(const libMesh::Elem &elem)
		{
			compute_affine_transformation(elem, A_inv_ , A_inv_m_b_);
		}

		void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override;
		void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

		AffineTransform3(const libMesh::DenseMatrix<libMesh::Real> &A_inv,
						 const libMesh::DenseVector<libMesh::Real> A_inv_m_b)
		: A_inv_(A_inv), A_inv_m_b_(A_inv_m_b)
		{}

		AffineTransform3()
		{}

		libMesh::DenseMatrix<libMesh::Real> &A_inv()
		{
			return A_inv_;
		}

		libMesh::DenseVector<libMesh::Real> &A_inv_m_b()
		{
			return A_inv_m_b_;
		}

	private:
		libMesh::DenseMatrix<libMesh::Real> A_inv_;
		libMesh::DenseVector<libMesh::Real> A_inv_m_b_;

		static void compute_affine_transformation(const libMesh::Elem &elem,
												  libMesh::DenseMatrix<libMesh::Real> &A_inv,
												  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
	};

	class SideAffineTransform3 : public Transform {
	public:
		inline SideAffineTransform3(const libMesh::Elem &elem, const int side)
		: a_trafo_()
		{
			compute_affine_transformation(elem, side, a_trafo_.A_inv(), a_trafo_.A_inv_m_b());
		}

		inline void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override
		{
			a_trafo_.transform_to_reference(world, ref);
			assert( std::abs(ref(2)) < 1e-8 );
		}

		void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

	private:
		AffineTransform3 a_trafo_;

		static void compute_affine_transformation(const libMesh::Elem &elem,
												  const int side,
												  libMesh::DenseMatrix<libMesh::Real> &A_inv,
												  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
	};

    class SideAffineTransform2 : public Transform {
    public:
        inline SideAffineTransform2(const libMesh::Elem &elem, const int side)
        : a_trafo_()
        {
            compute_affine_transformation(elem, side, a_trafo_.A_inv(), a_trafo_.A_inv_m_b());
        }

        inline void transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const override
        {
            a_trafo_.transform_to_reference(world, ref);

			//reference segment is (-1, 1)
			ref(0) *= 2.;
			ref(0) -= 1.;

            assert( std::abs(ref(1)) < 1e-8 );
            assert( std::abs(ref(2)) < 1e-8 );
        }

        void apply(const libMesh::Point &ref, libMesh::Point &world) const override { assert(false && "implement me"); }

    private:
        AffineTransform2 a_trafo_;

        static void compute_affine_transformation(const libMesh::Elem &elem,
                                                  const int side,
                                                  libMesh::DenseMatrix<libMesh::Real> &A_inv,
                                                  libMesh::DenseVector<libMesh::Real> &A_inv_m_b);
    };

	class QMortar : public libMesh::QBase {
	public:
		void resize(const int n_points)
		{
			this->get_points().resize(n_points);
			this->get_weights().resize(n_points);
		}

		QMortar(const unsigned int dim, const libMesh::Order order = libMesh::INVALID_ORDER)
		: libMesh::QBase(dim, order)
		{}

		virtual libMesh::QuadratureType type() const override
		{
			return libMesh::QGAUSS;
		}

		void init_1D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override
		{
			// assert(false);
		}

		void init_2D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override
		{
			// assert(false);
		}

		void init_3D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override
		{
			// assert(false);
		}
	};

	void print(const libMesh::QBase &ir, std::ostream &os = std::cout);
	double sum_of_weights(const libMesh::QBase &ir);
	double sum(const libMesh::DenseMatrix<libMesh::Real> &mat);

	void make_composite_quadrature_2D_non_affine(const libMesh::DenseMatrix<libMesh::Real> &polygon, const double weight, const int order, QMortar &c_ir);
	void make_composite_quadrature_2D_from_tri_mesh(const std::vector<int> &tri, const libMesh::DenseMatrix<libMesh::Real> &points,  const double weight, const int order, QMortar &c_ir);
	void make_composite_quadrature_2D(const libMesh::DenseMatrix<libMesh::Real> &polygon, const double weight, const int order, QMortar &c_ir);
	void make_composite_quadrature_3D(const Polyhedron &polyhedron, const double weight, const int order, QMortar &c_ir);
	void make_composite_quadrature_on_surf_2D(const libMesh::DenseMatrix<libMesh::Real> &line, const double weight, const int order, QMortar &c_ir);
	void make_composite_quadrature_on_surf_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon, const double weight, const int order, QMortar &c_ir);

	void transform_to_reference(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir);
	void transform_to_reference_surf(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir);
	double compute_volume(const Polyhedron &poly);

	void mortar_assemble(const libMesh::FEBase &trial_fe,
						 const libMesh::FEBase &test_fe,
						 libMesh::DenseMatrix<libMesh::Real> &elmat);

	void mortar_assemble(const libMesh::FEVectorBase &trial_fe,
						 const libMesh::FEVectorBase &test_fe,
						 libMesh::DenseMatrix<libMesh::Real> &elmat);

	void mortar_normal_and_gap_assemble(const libMesh::FEBase &test_fe,
										const libMesh::DenseVector<libMesh::Real> &surf_normal,
										const libMesh::DenseVector<libMesh::Real> &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap);

	void mortar_normal_and_gap_assemble(const uint dim,
										const libMesh::FEBase &test_fe,
										const libMesh::Point &surf_normal,
										const libMesh::Point &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap);

	void mortar_normal_and_gap_assemble(
										const libMesh::FEVectorBase &test_fe,
										const libMesh::DenseVector<libMesh::Real> &surf_normal,
										const libMesh::DenseVector<libMesh::Real> &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap);


	void mortar_normal_and_gap_assemble(const uint dim,
										const libMesh::FEVectorBase &test_fe,
										const libMesh::Point &surf_normal,
										const libMesh::Point &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap);

	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const libMesh::FEBase &test_fe,
											   const libMesh::DenseVector<libMesh::Real> &surf_normal,
											   const libMesh::DenseVector<libMesh::Real> &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap);

	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const uint dim,
											   const libMesh::FEBase &test_fe,
											   const libMesh::Point &surf_normal,
											   const libMesh::Point &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap);

	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const libMesh::FEVectorBase  &test_fe,
											   const libMesh::DenseVector<libMesh::Real> &surf_normal,
											   const libMesh::DenseVector<libMesh::Real> &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap);

	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const uint dim,
											   const libMesh::FEVectorBase  &test_fe,
											   const libMesh::Point &surf_normal,
											   const libMesh::Point &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap);

	// bool mortar_assemble(LibMeshFESpaceBase &src, LibMeshFESpaceBase &dest, std::shared_ptr<libMesh::SparseMatrix<libMesh::Real> > &B);

	// bool transfer(LibMeshFESpaceBase &src, libMesh::DenseVector<libMesh::Real> &src_fun, LibMeshFESpaceBase &dest, libMesh::DenseVector<libMesh::Real> &dest_fun);

	void make_polyhedron(const libMesh::Elem &e, Polyhedron &polyhedron);
	void make_polyline(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polyline);
	void make_polygon(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon);
	void make_polygon_3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon);

	bool intersect_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1, const libMesh::DenseMatrix<libMesh::Real> &poly2, libMesh::DenseMatrix<libMesh::Real> &intersection);
	bool intersect_3D(const libMesh::Elem &el1, const libMesh::Elem &el2, Polyhedron &intersection);
	bool intersect_3D(const Polyhedron &poly1, const Polyhedron &poly2, Polyhedron &intersection);

	bool project_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1,
					const libMesh::DenseMatrix<libMesh::Real> &poly2,
					libMesh::DenseMatrix<libMesh::Real> &projection_1,
					libMesh::DenseMatrix<libMesh::Real> &projection_2);

	bool project_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon_1,
					const libMesh::DenseMatrix<libMesh::Real> &polygon_2,
					libMesh::DenseMatrix<libMesh::Real> &projection_1,
					libMesh::DenseMatrix<libMesh::Real> &projection_2);

	bool biorthgonal_weights(const int type, libMesh::Real &w_ii, libMesh::Real &w_ij);

	void mortar_assemble_biorth(
								const libMesh::FEBase &trial_fe,
								const libMesh::FEBase &test_fe,
								const int type,
								libMesh::DenseMatrix<libMesh::Real> &elmat);

	void mortar_assemble_biorth(
								const libMesh::FEVectorBase &trial_fe,
								const libMesh::FEVectorBase &test_fe,
								const int type,
								libMesh::DenseMatrix<libMesh::Real> &elmat);


	void mortar_assemble_biorth(
								const int dim,
								const libMesh::FEBase &trial_fe,
								const libMesh::FEBase &test_fe,
								const int type,
								const libMesh::DenseVector<libMesh::Real> &indicator,
								libMesh::DenseMatrix<libMesh::Real> &elmat);

	void mortar_assemble_biorth(
								const int dim,
								const libMesh::FEVectorBase &trial_fe,
								const libMesh::FEVectorBase &test_fe,
								const int type,
								const libMesh::DenseVector<libMesh::Real> &indicator,
								libMesh::DenseMatrix<libMesh::Real> &elmat);

	void mortar_assemble_weights(const libMesh::FEVectorBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights);
	void mortar_assemble_weights(const libMesh::FEBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights);

	void mortar_assemble_weighted_biorth(
										 const libMesh::FEBase &trial_fe,
										 const libMesh::FEBase &test_fe,
										 const libMesh::DenseMatrix<libMesh::Real> &weights,
										 libMesh::DenseMatrix<libMesh::Real> &elmat);

	void mortar_assemble_weighted_biorth(
										 const libMesh::FEVectorBase &trial_fe,
										 const libMesh::FEVectorBase &test_fe,
										 const libMesh::DenseMatrix<libMesh::Real> &weights,
										 libMesh::DenseMatrix<libMesh::Real> &elmat);

	void mortar_normal_and_gap_assemble_weighted_biorth(
														const libMesh::FEVectorBase &test_fe,
														const int dim,
														const libMesh::Point &surf_normal,
														const libMesh::Point &plane_normal,
														const libMesh::Real &plane_offset,
														const libMesh::DenseMatrix<libMesh::Real> &weights,
														libMesh::DenseMatrix<libMesh::Real> &normals,
														libMesh::DenseVector<libMesh::Real> &gap);

	void mortar_normal_and_gap_assemble_weighted_biorth(
														const libMesh::FEBase &test_fe,
														const int dim,
														const libMesh::Point &surf_normal,
														const libMesh::Point &plane_normal,
														const libMesh::Real &plane_offset,
														const libMesh::DenseMatrix<libMesh::Real> &weights,
														libMesh::DenseMatrix<libMesh::Real> &normals,
														libMesh::DenseVector<libMesh::Real> &gap,
														const bool visdebug = false);
}

#endif //UTOPIA_MORTAR_ASSEMBLE_HPP
