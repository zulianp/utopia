
#include <libmesh/fe.h>
#include "utopia_triangulate.hpp"
#include "MortarAssemble.hpp"
#include "utopia_Polygon.hpp"
#include "utopia_intersector.hpp"

#include <memory>
#include <assert.h>
#include <algorithm>

namespace utopia {
	
	int order_for_l2_integral(const int dim,
							  const libMesh::Elem &master_el,
							  const int master_order,
							  const libMesh::Elem &slave_el,
							  const int slave_order)
	{
		bool slave_has_affine_map = slave_el.has_affine_map();
		
		int order = 0;
		if(dim == 2) {
			order = master_order * (is_quad(master_el.type())? 2 : 1 ) +
			slave_order  * (is_quad(slave_el.type()) ? 2 : 1 ) * (slave_has_affine_map? 1 : 2);
		} else if(dim == 3) {
			order = master_order * ( is_hex(master_el.type())? 3 : 1 ) +
			slave_order  * ( is_hex(slave_el.type())? 3  : 1 ) * (slave_has_affine_map? 1 : 2);
		} else {
			assert(false && "not supported yet for dim != 2 or dim != 3");
		}
		
		return order;
	}
	
	void Transform2::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const
	{
		ref = libMesh::FE<2, libMesh::LAGRANGE>::inverse_map(&elem_, world, 1e-10);
		assert( (libMesh::FE<2, libMesh::LAGRANGE>::on_reference_element(ref, elem_.type(), 1e-3)) );
		assert( (libMesh::FE<2, libMesh::LAGRANGE>::map(&elem_, ref).absolute_fuzzy_equals(world, 1e-8)) );
	}
	
	void AffineTransform2::compute_affine_transformation(const libMesh::Elem &elem, libMesh::DenseMatrix<libMesh::Real> &A_inv, libMesh::DenseVector<libMesh::Real> &A_inv_m_b)
	{
		libMesh::Point p0, p1, p2;
		
		switch(elem.type())
		{
			case libMesh::TRI3:
			case libMesh::TRI6:
			case libMesh::QUAD4:
			case libMesh::QUAD8:
			case libMesh::QUAD9:
			case libMesh::TRISHELL3:
			case libMesh::QUADSHELL4:
			// case libMesh::QUADSHELL8:
			{
				
				A_inv.resize(2,2);
				A_inv_m_b.resize(2);
				
				libMesh::DenseMatrix<libMesh::Real> A;
				
				A.resize(2,2);
				
				std::vector<const libMesh::Node *> elem_nodes;
				
				std::vector<libMesh::Point> reference_points;
				
				{
					libMesh::Point ref_p0(0.0, 0.0, 0.0);
					
					libMesh::Point ref_p1(1.0, 0.0, 0.0);
					
					libMesh::Point ref_p2(0.0, 1.0, 0.0);
					
					
					p0 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, ref_p0);
					
					p1 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, ref_p1);
					
					p2 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, ref_p2);
					
					
					
					A(0,0) =  (p1(0)-p0(0));
					A(0,1) =  (p2(0)-p0(0));
					A(1,0) =  (p1(1)-p0(1));
					A(1,1) =  (p2(1)-p0(1));
					
					
					libMesh::Real det =   A(0,0) * A(1,1) - A(0,1) * A(1,0);
					
					A_inv(0,0) = 1./det * A(1,1);
					A_inv(1,1) = 1./det * A(0,0);
					A_inv(0,1) = -1./det * A(0,1);
					A_inv(1,0) = -1./det * A(1,0);
					
					A_inv_m_b(0) = -1.0 * A_inv(0,0) * p0(0) - A_inv(0,1) * p0(1);
					A_inv_m_b(1) = -1.0 * A_inv(1,0) * p0(0) - A_inv(1,1) * p0(1);
					
				}
				
				break;
			}
			default:
			{
				assert(false && "implement me");
				break;
			}
		}
	}
	
	
	
	void AffineTransform2::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const
	{
		
		ref(0)= A_inv_m_b_(0) + A_inv_(0,0) * world(0) + A_inv_(0,1) * world(1);
		ref(1)= A_inv_m_b_(1) + A_inv_(1,0) * world(0) + A_inv_(1,1) * world(1);
		ref(2)= 0.0;
		
	}
	
	
	void Transform3::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const
	{
		ref = libMesh::FE<3, libMesh::LAGRANGE>::inverse_map(&elem_, world);
		// assert( (libMesh::FE<3, libMesh::LAGRANGE>::on_reference_element(ref, elem_.type(), 1e-6)) );
		assert( (libMesh::FE<3, libMesh::LAGRANGE>::map(&elem_, ref).absolute_fuzzy_equals(world, 1e-8)) );
	}
	
	
	
	void AffineTransform3::compute_affine_transformation(const libMesh::Elem &elem, libMesh::DenseMatrix<libMesh::Real> &A_inv, libMesh::DenseVector<libMesh::Real> &A_inv_m_b){
		
		
		libMesh::Point p0, p1, p2, p3;
		
		libMesh::Point ref_p0(0.0, 0.0, 0.0);
		libMesh::Point ref_p1(1.0, 0.0, 0.0);
		libMesh::Point ref_p2(0.0, 1.0, 0.0);
		libMesh::Point ref_p3(0.0, 0.0, 1.0);
		
		A_inv.resize(3,3);
		A_inv_m_b.resize(3);
		
		libMesh::DenseMatrix<libMesh::Real> A;
		
		A.resize(3,3);
		
		std::vector<const libMesh::Node *> elem_nodes;
		
		std::vector<libMesh::Point> reference_points;
		
		p0 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p0);
		p1 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p1);
		p2 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p2);
		p3 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, ref_p3);
		
		A(0,0) =  p1(0) - p0(0);
		A(0,1) =  p2(0) - p0(0);
		A(0,2) =  p3(0) - p0(0);
		A(1,0) =  p1(1) - p0(1);
		A(1,1) =  p2(1) - p0(1);
		A(1,2) =  p3(1) - p0(1);
		A(2,0) =  p1(2) - p0(2);
		A(2,1) =  p2(2) - p0(2);
		A(2,2) =  p3(2) - p0(2);
		
		libMesh::Real det =  Intersector::det_3(&A.get_values()[0]);
		Intersector::inverse_3(&A.get_values()[0], det, &A_inv.get_values()[0]);
		
		A_inv_m_b(0) = -1.0 * A_inv(0,0) * p0(0) - A_inv(0,1) * p0(1) - 1.0 * A_inv(0,2) * p0(2);
		A_inv_m_b(1) = -1.0 * A_inv(1,0) * p0(0) - A_inv(1,1) * p0(1) - 1.0 * A_inv(1,2) * p0(2);
		A_inv_m_b(2) = -1.0 * A_inv(2,0) * p0(0) - A_inv(2,1) * p0(1) - 1.0 * A_inv(2,2) * p0(2);
	}
	
	void SideAffineTransform2::compute_affine_transformation(const libMesh::Elem &elem,
															 const int side, libMesh::DenseMatrix<libMesh::Real> &A_inv,
															 libMesh::DenseVector<libMesh::Real> &A_inv_m_b)
	{
		//ref element -1, 1
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
		
		A_inv.resize(2,2);
		A_inv_m_b.resize(2);
		
		libMesh::DenseMatrix<libMesh::Real> A;
		A.resize(2, 2);
		
		A(0,0) = u(0);
		A(0,1) = n(0);
		A(1,0) = u(1);
		A(1,1) = n(1);
		
		const libMesh::Real det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
		
		A_inv(0,0) =  1./det * A(1, 1);
		A_inv(1,1) =  1./det * A(0, 0);
		A_inv(0,1) = -1./det * A(0, 1);
		A_inv(1,0) = -1./det * A(1, 0);
		
		A_inv_m_b(0) = -1.0 * A_inv(0, 0) * p0(0) - A_inv(0, 1) * p0(1);
		A_inv_m_b(1) = -1.0 * A_inv(1, 0) * p0(0) - A_inv(1, 1) * p0(1);		
	}
	
	void SideAffineTransform3::compute_affine_transformation(const libMesh::Elem &elem,  const int side, libMesh::DenseMatrix<libMesh::Real> &A_inv, libMesh::DenseVector<libMesh::Real> &A_inv_m_b)
	{
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
		
		A_inv.resize(3,3);
		A_inv_m_b.resize(3);
		
		libMesh::DenseMatrix<libMesh::Real> A;
		
		A.resize(3,3);
		
		A(0,0) = u(0);
		A(0,1) = v(0);
		A(0,2) = n(0);
		A(1,0) = u(1);
		A(1,1) = v(1);
		A(1,2) = n(1);
		A(2,0) = u(2);
		A(2,1) = v(2);
		A(2,2) = n(2);
		
		libMesh::Real det =  Intersector::det_3(&A.get_values()[0]);
		Intersector::inverse_3(&A.get_values()[0], det, &A_inv.get_values()[0]);
		
		A_inv_m_b(0) = -1.0 * A_inv(0,0) * p0(0) - A_inv(0,1) * p0(1) - 1.0 * A_inv(0,2) * p0(2);
		A_inv_m_b(1) = -1.0 * A_inv(1,0) * p0(0) - A_inv(1,1) * p0(1) - 1.0 * A_inv(1,2) * p0(2);
		A_inv_m_b(2) = -1.0 * A_inv(2,0) * p0(0) - A_inv(2,1) * p0(1) - 1.0 * A_inv(2,2) * p0(2);
	}
	
	void AffineTransform3::transform_to_reference(const libMesh::Point &world, libMesh::Point &ref) const
	{
		ref(0) = A_inv_m_b_(0) + A_inv_(0,0) * world(0) + A_inv_(0,1) * world(1) + A_inv_(0,2) * world(2);
		ref(1) = A_inv_m_b_(1) + A_inv_(1,0) * world(0) + A_inv_(1,1) * world(1) + A_inv_(1,2) * world(2);
		ref(2) = A_inv_m_b_(2) + A_inv_(2,0) * world(0) + A_inv_(2,1) * world(1) + A_inv_(2,2) * world(2);
	}
	
	void print(const libMesh::QBase &ir, std::ostream &os)
	{
		os << "points:\n[\n";
		for(int i = 0; i < ir.n_points(); ++i) {
			os << "\t" << ir.qp(i)(0) << ", " << ir.qp(i)(1) << ", " << ir.qp(i)(2) << "\n";
		}
		
		os << "]\n";
		
		os << "weights:\n[\n\t";
		for(int i = 0; i < ir.n_points(); ++i) {
			os << ir.w(i) << " ";
		}
		os << "\n]\n";
	}
	
	double sum_of_weights(const libMesh::QBase &ir)
	{
		double ret = 0;
		for(int i = 0; i < ir.n_points(); ++i) {
			ret += ir.w(i);
		}
		return ret;
	}
	
	template<class V, class T>
	static void add(const V &p, const T &alpha, const V &v, V &result)
	{
		assert(result.size() == p.size());
		
		for(int i = 0; i < result.size(); ++i) {
			result(i) = p(i) + alpha * v(i);
		}
	}
	
	void make_composite_quadrature_2D_non_affine(const libMesh::DenseMatrix<libMesh::Real> &polygon, const double weight, const int order, QMortar &c_ir)
	{
		std::vector<int> tri;
		//triangulate_polygon(polygon.m(), &polygon.get_values()[0], tri);
		//make_composite_quadrature_2D_from_tri_mesh(tri, polygon, weight, order, c_ir);
		assert(false);
	}
	
	
	void make_composite_quadrature_2D_from_tri_mesh(const std::vector<int> &tri, const libMesh::DenseMatrix<libMesh::Real> &points,  const double weight, const int order, QMortar &c_ir)
	{
		using namespace libMesh;
		//FIXME only for sligthlty deformed elements
		libMesh::QGauss ir(2, libMesh::Order(order));
		ir.init(libMesh::TRI6);
		
		double triangle[3 * 2] = { 0., 0., 0., 0., 0., 0. };
		
		const int n_triangles   = tri.size() / 3;
		const int n_quad_points = n_triangles * ir.n_points();
		c_ir.resize(n_quad_points);
		
		
		libMesh::DenseVector<libMesh::Real> o, u, v, p;
		double relative_weight = 0;
		int quad_index = 0;
		for(int i = 0; i < n_triangles; ++i) {
			const int i3 = i * 3;
			int v1 = tri[i3];
			int v2 = tri[i3 + 1];
			int v3 = tri[i3 + 2];
			
			get_row(v1, points, o);
			get_row(v2, points, u);
			get_row(v3, points, v);
			
			triangle[0] = o(0);
			triangle[1] = o(1);
			
			triangle[2] = u(0);
			triangle[3] = u(1);
			
			triangle[4] = v(0);
			triangle[5] = v(1);
			
			//plot_polygon(2, 3, triangle);
			
			u -= o;
			v -= o;
			
			const double scale = fabs(Intersector::polygon_area_2(3, triangle)) / ( weight );
			relative_weight += scale;
			
			for(int k = 0; k < ir.n_points(); ++k, ++quad_index) {
				auto &qp    = ir.get_points()[k];
				auto &c_qp 	= c_ir.get_points()[quad_index];
				
				p = o;
				add(p, qp(0), u, p);
				add(p, qp(1), v, p);
				
				c_qp(0) = p(0);
				c_qp(1) = p(1);
				c_qp(2) = 0.0;
				
				c_ir.get_weights()[quad_index] = ir.w(k) * scale;
			}
		}
		
		assert(relative_weight <= 1.0001);
		assert(quad_index == n_quad_points);
		// plot_quad_points(2, c_ir.get_points(), "qp");
	}
	
	void make_composite_quadrature_2D(const libMesh::DenseMatrix<libMesh::Real> &polygon, const double weight, const int order, QMortar &c_ir)
	{
		libMesh::QGauss ir(2, libMesh::Order(order));
		ir.init(libMesh::TRI3);
		
		// std::cout << "ref quad " << std::endl;
		// ir.print_info();
		// std::cout << "```````````" << std::endl;
		
		const int n_triangles     = polygon.m() - 2;
		const int n_quad_points   = n_triangles * ir.n_points();
		
		assert(fabs(sum_of_weights(ir) - 0.5) < 1e-8);
		
		double triangle[3*2] = { polygon.get_values()[0], polygon.get_values()[1], 0., 0., 0., 0. };
		
		libMesh::DenseVector<libMesh::Real> o, u, v, p;
		get_row(0, polygon, o);
		
		c_ir.resize(n_quad_points);
		
		double relative_weight = 0;
		
		int quad_index = 0;
		for(int i = 2; i < polygon.m(); ++i) {
			get_row(i-1, polygon, u);
			get_row(i,   polygon, v);
			
			triangle[2] = u(0);
			triangle[3] = u(1);
			
			triangle[4] = v(0);
			triangle[5] = v(1);
			
			
			u -= o;
			v -= o;
			
			const double scale = fabs(Intersector::polygon_area_2(3, triangle)) / ( weight );
			relative_weight += scale;
			
			for(int k = 0; k < ir.n_points(); ++k, ++quad_index) {
				auto &qp    = ir.get_points()[k];
				auto &c_qp 	= c_ir.get_points()[quad_index];
				
				p = o;
				add(p, qp(0), u, p);
				add(p, qp(1), v, p);
				
				c_qp(0) = p(0);
				c_qp(1) = p(1);
				c_qp(2) = 0.0;
				
				c_ir.get_weights()[quad_index] = ir.w(k) * scale;
			}
		}
		
		assert(relative_weight <= 1.0001);
		assert(quad_index == n_quad_points);
	}
	
	void make_composite_quadrature_3D(const Polyhedron &polyhedron, const double weight, const int order, QMortar &c_ir)
	{
		using std::min;
		
		
		libMesh::QGauss ir(3, libMesh::Order(order));
		ir.init(libMesh::TET4);
		
		const int n_sub_elements      = Intersector::n_volume_elements(polyhedron);
		const int total_n_quad_points = n_sub_elements * ir.n_points();
		
		std::vector<double> ref_quad_points (ir.n_points() * 3);
		std::vector<double> ref_quad_weights(ir.n_points());
		
		for(int i = 0; i < ir.n_points(); ++i) {
			const int offset = i * 3;
			
			ref_quad_points[offset 	  ] = ir.qp(i)(0);
			ref_quad_points[offset + 1] = ir.qp(i)(1);
			ref_quad_points[offset + 2] = ir.qp(i)(2);
			
			ref_quad_weights[i] 		= ir.w(i);
		}
		
		c_ir.resize(total_n_quad_points);
		
		//global quadrature points
		double quad_points     [MAX_QUAD_POINTS * 3];
		double quad_weights    [MAX_QUAD_POINTS];
		
		double barycenter_p[3];
		Intersector::row_average( polyhedron.n_nodes, polyhedron.n_dims, polyhedron.points, barycenter_p);
		
		const uint max_n_sub_els = Intersector::max_n_elements_from_facets(polyhedron);
		const uint max_n_sub_inc = std::max(1u, MAX_QUAD_POINTS / (max_n_sub_els * ir.n_points()));
		
		const double scale = 1.0 / ( weight );
		
		int utopia_fe_quad_index = 0;
		
		if(n_sub_elements == 1) {
			Intersector::tetrahedron_transform(polyhedron.points, total_n_quad_points, &ref_quad_points[0], quad_points);
			
			const double w = fabs(Intersector::m_tetrahedron_volume(polyhedron.points) * scale);
			for(int i = 0; i < total_n_quad_points; ++i, ++utopia_fe_quad_index) {
				const int offset = i * 3;
				c_ir.get_points()[utopia_fe_quad_index](0) = quad_points[offset 	];
				c_ir.get_points()[utopia_fe_quad_index](1) = quad_points[offset + 1];
				c_ir.get_points()[utopia_fe_quad_index](2) = quad_points[offset + 2];
				
				c_ir.get_weights()[utopia_fe_quad_index] = ref_quad_weights[i] * w;
			}
			
		} else {
			
			for(int begin_k = 0; begin_k < polyhedron.n_elements;) {
				const int end_k = min(begin_k + max_n_sub_inc, static_cast<uint>(polyhedron.n_elements));
				assert(end_k > begin_k && "end_k > begin_k");
				
				const int n_quad_points =
				Intersector::make_quadrature_points_from_polyhedron_in_range_around_point(
																					 polyhedron,
																					 begin_k,
																					 end_k,
																					 scale,
																					 ir.n_points(),
																					 &ref_quad_points [0],
																					 &ref_quad_weights[0],
																					 barycenter_p,
																					 quad_points,
																					 quad_weights
																					 );
				
				for(int i = 0; i < n_quad_points; ++i, ++utopia_fe_quad_index) {
					const int offset = i * 3;
					c_ir.get_points()[utopia_fe_quad_index](0) = quad_points[offset    ];
					c_ir.get_points()[utopia_fe_quad_index](1) = quad_points[offset + 1];
					c_ir.get_points()[utopia_fe_quad_index](2) = quad_points[offset + 2];
					
					c_ir.get_weights()[utopia_fe_quad_index] = quad_weights[i];
				}
				
				begin_k = end_k;
			}
		}
	}
	
	void make_composite_quadrature_on_surf_2D(const libMesh::DenseMatrix<libMesh::Real> &line, const double weight, const int order, QMortar &c_ir)
	{
		using namespace libMesh;
		QGauss ir(1, libMesh::Order(order));
		ir.init(libMesh::EDGE2);
		
		c_ir.resize(ir.n_points());
		
		Point o(2), v(2);
		
		o(0) = line(0, 0);
		o(1) = line(0, 1);
		
		v(0) = line(1, 0) - line(0, 0);
		v(1) = line(1, 1) - line(0, 1);

		const Real length = v.norm();
		
		for(int k = 0; k < ir.n_points(); ++k) {
			const Real w = (1 + ir.get_points()[k](0)) * 0.5;
			assert(w >= -1e-16);
			assert(w <= 1 + 1e-16);
			
			c_ir.get_points()[k] = v;
			c_ir.get_points()[k] *= w;
			c_ir.get_points()[k] += o;
			
			c_ir.get_weights()[k] = ir.w(k) * length * weight * 0.5;
		}

	}
	
	void make_composite_quadrature_on_surf_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon, const double weight, const int order, QMortar &c_ir)
	{
		
		libMesh::QGauss ir(2, libMesh::Order(order));

		if(order <= 2) {
			ir.init(libMesh::TRI3);
		} else {
			ir.init(libMesh::TRI6);
		}
		
		// std::cout << "ref quad " << std::endl;
		// ir.print_info();
		// std::cout << "```````````" << std::endl;
		
		const int n_triangles     = polygon.m() - 2;
		const int n_quad_points   = n_triangles * ir.n_points();
		
		assert(fabs(sum_of_weights(ir) - 0.5) < 1e-8);
		
		double triangle[3 * 3] = {  polygon.get_values()[0], polygon.get_values()[1],  polygon.get_values()[2],
			0., 0., 0.,
			0., 0., 0.
		};
		
		libMesh::DenseVector<libMesh::Real> o, u, v, p;
		get_row(0, polygon, o);
		
		c_ir.resize(n_quad_points);
		
		double relative_weight = 0;
		
		int quad_index = 0;
		for(int i = 2; i < polygon.m(); ++i) {
			get_row(i-1, polygon, u);
			get_row(i,   polygon, v);
			
			triangle[3] = u(0);
			triangle[4] = u(1);
			triangle[5] = u(2);
			
			triangle[6] = v(0);
			triangle[7] = v(1);
			triangle[8] = v(2);
			
			u -= o;
			v -= o;
			
			relative_weight = Intersector::polygon_area_3(3, triangle) * weight * 2.; 
			assert(relative_weight >= 0.);
			
			for(int k = 0; k < ir.n_points(); ++k, ++quad_index) {
				auto &qp    = ir.get_points()[k];
				auto &c_qp 	= c_ir.get_points()[quad_index];
				
				p = o;
				add(p, qp(0), u, p);
				add(p, qp(1), v, p);
				
				c_qp(0) = p(0);
				c_qp(1) = p(1);
				c_qp(2) = p(2);
				
				c_ir.get_weights()[quad_index] = ir.w(k) * relative_weight;
			}
		}
		
		assert(quad_index == n_quad_points);
	}
	
	void transform_to_reference(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir)
	{
		const int dim = global_ir.get_dim();
		libMesh::Point p;
		ref_ir.resize(global_ir.n_points());
		
		libMesh::DenseMatrix<libMesh::Real> A_inv;
		
		for(int i = 0; i < global_ir.n_points(); ++i) {
			p(0) = global_ir.qp(i)(0);
			p(1) = global_ir.qp(i)(1);
			
			if(dim > 2) {
				p(2) = global_ir.qp(i)(2);
			}
			
			
			trans.transform_to_reference(p, ref_ir.get_points()[i]);
			
			ref_ir.get_weights()[i] = global_ir.w(i);
			
			if(is_hex(type)) {
				//remove tet scaling and add hex scaling
				ref_ir.get_weights()[i] *= (8.*6.);
			} else if(is_quad(type)) {
				//remove triangle scaling and add quad scaling
				ref_ir.get_weights()[i] *= (4.*2.);
			}
		}
	}
	
	void transform_to_reference_surf(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir)
	{
		assert(is_valid_elem_type(type));
		
		const int dim = global_ir.get_dim();
		libMesh::Point p;
		ref_ir.resize(global_ir.n_points());
		
		
		libMesh::DenseMatrix<libMesh::Real> A_inv;
		
		for(int i = 0; i < global_ir.n_points(); ++i) {
			p(0) = global_ir.qp(i)(0);
			p(1) = global_ir.qp(i)(1);
			
			if(dim > 2) {
				p(2) = global_ir.qp(i)(2);
			}
			
			trans.transform_to_reference(p, ref_ir.get_points()[i]);
			
			ref_ir.get_weights()[i] = global_ir.w(i);
			
			
			
			if(is_tri(type)) {
				ref_ir.get_weights()[i] *= 2.;
			} else if(is_quad(type)) {
				ref_ir.get_weights()[i] *= 2.;
			} else if(is_hex(type)) {
				ref_ir.get_weights()[i] *= 1.;
			} else if(is_tet(type)) {
				ref_ir.get_weights()[i] *= 0.5;
			} else {
				assert(false && "add special case");
			}
		}
	}
	
	template<class Left, class Right>
	libMesh::Real contract(const Left &left, const Right &right)
	{
		return left.contract(right);
	}
	
	libMesh::Real contract(const libMesh::Real &left, const libMesh::Real &right)
	{
		return left * right;
	}
	
	template<class FE>
	void mortar_assemble_aux(
							 const FE &trial_fe,
							 const FE &test_fe,
							 libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		if(elmat.m() != test_fe.get_phi().size() ||  elmat.n() != trial_fe.get_phi().size()) {
			elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
			elmat.zero();
		}
		
		const auto &trial = trial_fe.get_phi();
		const auto &test  = test_fe.get_phi();
		const auto &JxW   = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_trial = trial.size();
		const uint n_qp    = test[0].size();
		
		assert(test[0].size() == trial[0].size());
		
		for(uint i = 0; i < n_test; ++i) {
			for(uint j = 0; j < n_trial; ++j) {
				for(uint qp = 0; qp < n_qp; ++qp) {
					elmat(i, j) += contract(test[i][qp], trial[j][qp]) * JxW[qp];
				}
			}
		}
	}
	
	
	// template<class FE>
	// void mortar_assemble_aux_reverse(
	//                          const FE &test_fe,
	//                          const FE &trial_fe,
	//                          libMesh::DenseMatrix<libMesh::Real> &elmat)
	// {
	//     if(elmat.m() != test_fe.get_phi().size() ||  elmat.n() != trial_fe.get_phi().size()) {
	//         elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
	//         elmat.zero();
	//     }
	
	//     const auto &trial = trial_fe.get_phi();
	//     const auto &test  = test_fe.get_phi();
	//     const auto &JxW   = test_fe.get_JxW();
	
	//     const uint n_test  = test.size();
	//     const uint n_trial = trial.size();
	//     const uint n_qp    = test[0].size();
	
	//     for(uint i = 0; i < n_test; ++i) {
	//         for(uint j = 0; j < n_trial; ++j) {
	//             for(uint qp = 0; qp < n_qp; ++qp) {
	//                 elmat(i, j) += contract(test[i][qp], trial[j][qp]) * JxW[qp];
	//             }
	//         }
	//     }
	// }
	
	
	
	libMesh::Real len(const libMesh::Real val)
	{
		return std::abs(val);
	}
	
	template<class Vec>
	libMesh::Real len(const Vec &val)
	{
		return val.size();
	}
	
	static inline bool is_vec(const libMesh::Real val)
	{
		return false;
	}
	
	template<class Vec>
	bool is_vec(const Vec &val)
	{
		return true;
	}
	
	void make_tp(const int i,  libMesh::Real &val) {}
	
	template<class Vec>
	void make_tp(const int i, Vec &val) {
		libMesh::Real s =  val(i);
		val.zero();
		val(i) = s;
	}
	
	
	template<class FE>
	void mortar_assemble_biorth_aux(
									const FE &trial_fe,
									const FE &test_fe,
									const libMesh::Real &wii,
									const libMesh::Real &wij,
									libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		if(elmat.m() != test_fe.get_phi().size() ||  elmat.n() != trial_fe.get_phi().size()) {
			elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
			elmat.zero();
		}
		
		const auto &trial = trial_fe.get_phi();
		const auto &test  = test_fe.get_phi();
		const auto &JxW   = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_trial = trial.size();
		const uint n_qp    = test[0].size();
		
		for(uint qp = 0; qp < n_qp; ++qp) {
			for(uint i = 0; i < n_test; ++i) {
				
				auto biorth_test = ((0 == i) ? wii : wij) * test[0][qp];
				
				for(uint k = 1; k < n_test; ++k) {
					biorth_test += ((k == i) ? wii : wij) * test[k][qp];
				}
				
				for(uint j = 0; j < n_trial; ++j) {
					elmat(i, j) += contract(biorth_test, trial[j][qp]) * JxW[qp];
				}
			}
		}
	}
	
	void mortar_assemble(
						 const libMesh::FEBase &trial_fe,
						 const libMesh::FEBase &test_fe,
						 libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		mortar_assemble_aux(trial_fe, test_fe, elmat);
	}
	
	void mortar_assemble(const libMesh::FEVectorBase &trial_fe,
						 const libMesh::FEVectorBase &test_fe,
						 libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		mortar_assemble_aux(trial_fe, test_fe, elmat);
	}
	
	void mortar_assemble_biorth(
								const libMesh::FEBase &trial_fe,
								const libMesh::FEBase &test_fe,
								const int type,
								libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		
		
		libMesh::Real w_ii, w_ij;
		biorthgonal_weights(type, w_ii, w_ij);
		mortar_assemble_biorth_aux(trial_fe, test_fe, w_ii, w_ij, elmat);
	}
	
	void mortar_assemble_biorth(
								const libMesh::FEVectorBase &trial_fe,
								const libMesh::FEVectorBase &test_fe,
								const int type,
								libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		libMesh::Real w_ii, w_ij;
		biorthgonal_weights(type, w_ii, w_ij);
		mortar_assemble_biorth_aux(trial_fe, test_fe, w_ii, w_ij, elmat);
	}
	
	
	template<class FE>
	void mortar_assemble_biorth_aux(
									const int dim,
									const FE &trial_fe,
									const FE &test_fe,
									const libMesh::Real &wii,
									const libMesh::Real &wij,
									const libMesh::DenseVector<libMesh::Real> &indicator,
									libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		if(elmat.m() != test_fe.get_phi().size() ||  elmat.n() != trial_fe.get_phi().size()) {
			elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
			elmat.zero();
		}
		
		const auto &trial = trial_fe.get_phi();
		const auto &test  = test_fe.get_phi();
		const auto &JxW   = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_trial = trial.size();
		const uint n_qp    = test[0].size();
		
//		bool v = is_vec(test[0][0]);
		
		for(uint qp = 0; qp < n_qp; ++qp) {
			for(uint i = 0; i < n_test; ++i) {
				
				auto biorth_test = ((0 == i) ? wii : wij) * test[0][qp];
				
				for(uint k = 1; k < n_test; ++k) {
					biorth_test += ((k == i) ? wii : wij) * test[k][qp];
				}
				
				for(int k = 0; k < dim; ++k) {
					if(i % dim == k) {
						make_tp(k, biorth_test);
						break;
					}
				}
				
				// if(indicator(i) > 0)
				// 	std::cout <<  biorth_test << std::endl;
				
				for(uint j = 0; j < n_trial; ++j) {
					elmat(i, j) +=  indicator(i) * contract(biorth_test, trial[j][qp]) * JxW[qp];
				}
			}
		}
	}
	
	
	void mortar_assemble_biorth(
								const int dim,
								const libMesh::FEBase &trial_fe,
								const libMesh::FEBase &test_fe,
								const int type,
								const libMesh::DenseVector<libMesh::Real> &indicator,
								libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		
		
		libMesh::Real w_ii, w_ij;
		biorthgonal_weights(type, w_ii, w_ij);
		mortar_assemble_biorth_aux(dim, trial_fe, test_fe, w_ii, w_ij, indicator, elmat);
	}
	
	void mortar_assemble_biorth(
								const int dim,
								const libMesh::FEVectorBase &trial_fe,
								const libMesh::FEVectorBase &test_fe,
								const int type,
								const libMesh::DenseVector<libMesh::Real> &indicator,
								libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		libMesh::Real w_ii, w_ij;
		biorthgonal_weights(type, w_ii, w_ij);
		mortar_assemble_biorth_aux(dim, trial_fe, test_fe, w_ii, w_ij, indicator, elmat);
	}
	
	
	void mortar_normal_and_gap_assemble(const uint dim,
										const libMesh::FEBase &test_fe,
										const libMesh::Point &surf_normal,
										const libMesh::Point &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap)
	{
		using namespace libMesh;
		DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);
		
		for(uint i = 0; i < dim; ++i) {
			surf_normal_v(i)  = surf_normal(i);
			plane_normal_v(i) = plane_normal(i);
		}
		
		mortar_normal_and_gap_assemble(test_fe, surf_normal_v, plane_normal_v, plane_offset, normals, gap);
	}
	
	void mortar_normal_and_gap_assemble(
										const libMesh::FEBase &test_fe,
										const libMesh::DenseVector<libMesh::Real> &surf_normal,
										const libMesh::DenseVector<libMesh::Real> &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap)
	{
		using namespace libMesh;
		
		
		const uint dim = plane_normal.size();
		
		if(normals.m() != test_fe.get_phi().size() || dim != normals.n()) {
			normals.resize(test_fe.get_phi().size(), dim);
			normals.zero();
			gap.resize(test_fe.get_phi().size());
			gap.zero();
		}
		
		const auto &test   = test_fe.get_phi();
		// const auto &grad   = test_fe.get_dphi();
		const auto &point  = test_fe.get_xyz();
		const auto &JxW    = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_qp    = test[0].size();
		
		DenseVector<Real> p(dim);
		
		for(uint i = 0; i < n_test; ++i) {
			for(uint qp = 0; qp < n_qp; ++qp) {
				
				p(0) = point[qp](0);
				p(1) = point[qp](1);
				
				if(dim > 2) {
					p(2) = point[qp](2);
				}
				
				Real isect = 0;
				Intersector::intersect_ray_with_plane(dim, 1, &p.get_values()[0], &surf_normal.get_values()[0], &plane_normal.get_values()[0], plane_offset, &isect);
				
				
				// printf("g: %g (%g, %g)\n", isect, p(1), plane_offset);
				// assert(isect > 0);
				
				gap(i) += test[i][qp] * isect * JxW[qp];
				
				for(uint d = 0; d < dim; ++d) {
					normals(i, d) += test[i][qp] * surf_normal(d) * JxW[qp];
				}
			}
		}
	}
	
	void mortar_normal_and_gap_assemble(const uint dim,
										const libMesh::FEVectorBase &test_fe,
										const libMesh::Point &surf_normal,
										const libMesh::Point &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap)
	{
		using namespace libMesh;
		DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);
		
		for(uint i = 0; i < dim; ++i) {
			surf_normal_v(i)  = surf_normal(i);
			plane_normal_v(i) = plane_normal(i);
		}
		
		mortar_normal_and_gap_assemble(test_fe, surf_normal_v, plane_normal_v, plane_offset, normals, gap);
	}
	
	void mortar_normal_and_gap_assemble(
										const libMesh::FEVectorBase &test_fe,
										const libMesh::DenseVector<libMesh::Real> &surf_normal,
										const libMesh::DenseVector<libMesh::Real> &plane_normal,
										const libMesh::Real &plane_offset,
										libMesh::DenseMatrix<libMesh::Real> &normals,
										libMesh::DenseVector<libMesh::Real> &gap)
	{
		using namespace libMesh;
		
		
		const uint dim = plane_normal.size();
		
		if(normals.m() != test_fe.get_phi().size()/dim || dim != normals.n()) {
			normals.resize(test_fe.get_phi().size()/dim, dim);
			normals.zero();
			gap.resize(test_fe.get_phi().size());
			gap.zero();
		}
		
		const auto &test   = test_fe.get_phi();
		const auto &point  = test_fe.get_xyz();
		const auto &JxW    = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_qp    = test[0].size();
		
		DenseVector<Real> p(dim);
		DenseVector<Real> v(dim);
		
		for(uint qp = 0; qp < n_qp; ++qp) {
			
			p(0) = point[qp](0);
			p(1) = point[qp](1);
			
			if(dim > 2) {
				p(2) = point[qp](2);
			}
			
			Real isect = 0;
			Intersector::intersect_ray_with_plane(dim, 1, &p.get_values()[0], &surf_normal.get_values()[0], &plane_normal.get_values()[0], plane_offset, &isect);
			
			v = surf_normal;
			v *= isect;
			// quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);
			
			for(uint i = 0; i < n_test; ++i) {
				gap(i) += test[i][qp](0) * isect * JxW[qp];
				
				for(uint d = 0; d < dim; ++d) {
					normals.get_values()[i] += test[i][qp](d) * surf_normal(d) * JxW[qp];
				}
			}
		}
	}
	
	
	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const uint dim,
											   const libMesh::FEBase &test_fe,
											   const libMesh::Point &surf_normal,
											   const libMesh::Point &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap)
	{
		using namespace libMesh;
		DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);
		
		for(uint i = 0; i < dim; ++i) {
			surf_normal_v(i)  = surf_normal(i);
			plane_normal_v(i) = plane_normal(i);
		}
		
		mortar_normal_and_gap_assemble_biorth(type, test_fe, surf_normal_v, plane_normal_v, plane_offset, indicator, normals, gap);
	}
	
	
	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const libMesh::FEBase &test_fe,
											   const libMesh::DenseVector<libMesh::Real> &surf_normal,
											   const libMesh::DenseVector<libMesh::Real> &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap)
	{
		using namespace libMesh;
		
		
		const uint dim = plane_normal.size();
		
		libMesh::Real w_ii, w_ij;
		biorthgonal_weights(type, w_ii, w_ij);
		
		if(normals.m() != test_fe.get_phi().size() || dim != normals.n()) {
			normals.resize(test_fe.get_phi().size(), dim);
			normals.zero();
			gap.resize(test_fe.get_phi().size());
			gap.zero();
		}
		
		const auto &test   = test_fe.get_phi();
		// const auto &grad   = test_fe.get_dphi();
		const auto &point  = test_fe.get_xyz();
		const auto &JxW    = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_qp    = test[0].size();
		
		DenseVector<Real> p(dim);
		// DenseVector<Real> v(dim);
		
		for(uint i = 0; i < n_test; ++i) {
			for(uint qp = 0; qp < n_qp; ++qp) {
				
				p(0) = point[qp](0);
				p(1) = point[qp](1);
				
				if(dim > 2) {
					p(2) = point[qp](2);
				}
				
				Real isect = 0;
				Intersector::intersect_ray_with_plane(dim, 1, &p.get_values()[0], &surf_normal.get_values()[0], &plane_normal.get_values()[0], plane_offset, &isect);
				
				// v = surf_normal;
				// v *= isect;
				// quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);
				
				// printf("g: %g (%g, %g)\n", isect, p(1), plane_offset);
				// assert(isect > 0);
				
				auto biorth_test = ((0 == i) ? w_ii : w_ij) * test[0][qp];
				
				for(uint k = 1; k < n_test; ++k) {
					biorth_test += ((k == i) ? w_ii : w_ij) * test[k][qp];
				}
				
				for(int k = 0; k < dim; ++k) {
					if(i % dim == k) {
						make_tp(k, biorth_test);
						break;
					}
				}
				
				
				gap(i) += indicator(i) * biorth_test * isect * JxW[qp];
				
				for(uint d = 0; d < dim; ++d) {
					normals(i, d) += indicator(i) * biorth_test * surf_normal(d) * JxW[qp];
				}
			}
		}
	}
	
	
	
	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const libMesh::FEVectorBase &test_fe,
											   const libMesh::DenseVector<libMesh::Real> &surf_normal,
											   const libMesh::DenseVector<libMesh::Real> &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap)
	{
		libMesh::Real w_ii, w_ij;
		biorthgonal_weights(type, w_ii, w_ij);
		
		// std::cout << w_ii << " " << w_ij << std::endl;
		
		using namespace libMesh;
		
		
		const uint dim = plane_normal.size();
		
		if(normals.m() != test_fe.get_phi().size()/dim || dim != normals.n()) {
			normals.resize(test_fe.get_phi().size()/dim, dim);
			normals.zero();
			gap.resize(test_fe.get_phi().size());
			gap.zero();
		}
		
		const auto &test   = test_fe.get_phi();
		const auto &point  = test_fe.get_xyz();
		const auto &JxW    = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_qp    = test[0].size();
		
		DenseVector<Real> p(dim);
		DenseVector<Real> v(dim);
		
		for(uint qp = 0; qp < n_qp; ++qp) {
			
			p(0) = point[qp](0);
			p(1) = point[qp](1);
			
			if(dim > 2) {
				p(2) = point[qp](2);
			}
			
			Real isect = 0;
			Intersector::intersect_ray_with_plane(dim, 1, &p.get_values()[0], &surf_normal.get_values()[0], &plane_normal.get_values()[0], plane_offset, &isect);
			
			v = surf_normal;
			v *= isect;
			// quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);
			
			for(uint i = 0; i < n_test; ++i) {
				auto biorth_test = ((0 == i) ? w_ii : w_ij) * test[0][qp];
				
				for(uint k = 1; k < n_test; ++k) {
					biorth_test += ((k == i) ? w_ii : w_ij) * test[k][qp];
				}
				
				for(int k = 0; k < dim; ++k) {
					if(i % dim == k) {
						make_tp(k, biorth_test);
						break;
					}
				}
				
				gap(i) +=  indicator(i) * biorth_test(0) * isect * JxW[qp];
				
				for(uint d = 0; d < dim; ++d) {
					normals.get_values()[i] +=  indicator(i) * biorth_test(d) * surf_normal(d) * JxW[qp];
				}
			}
		}
	}
	
	void mortar_normal_and_gap_assemble_biorth(
											   const int type,
											   const uint dim,
											   const libMesh::FEVectorBase &test_fe,
											   const libMesh::Point &surf_normal,
											   const libMesh::Point &plane_normal,
											   const libMesh::Real &plane_offset,
											   const libMesh::DenseVector<libMesh::Real> &indicator,
											   libMesh::DenseMatrix<libMesh::Real> &normals,
											   libMesh::DenseVector<libMesh::Real> &gap)
	{
		using namespace libMesh;
		DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);
		
		for(uint i = 0; i < dim; ++i) {
			surf_normal_v(i)  = surf_normal(i);
			plane_normal_v(i) = plane_normal(i);
		}
		
		mortar_normal_and_gap_assemble_biorth(type, test_fe, surf_normal_v, plane_normal_v, plane_offset, indicator, normals, gap);
	}
	
	bool intersect_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1, const libMesh::DenseMatrix<libMesh::Real> &poly2, libMesh::DenseMatrix<libMesh::Real> &intersection)
	{
		double result_buffer[MAX_N_ISECT_POINTS * 2];
		int n_vertices_result;
		
		assert( Intersector::polygon_area_2(poly1.m(),  &poly1.get_values()[0]) > 0 );
		assert( Intersector::polygon_area_2(poly2.m(),  &poly2.get_values()[0]) > 0 );
		
		if(!
			// Intersector::
			intersect_convex_polygons(poly1.m(),
									  &poly1.get_values()[0],
									  poly2.m(),
									  &poly2.get_values()[0],
									  &n_vertices_result,
									  result_buffer,
									  DEFAULT_TOLLERANCE)) {
			return false;
		}
		
		assert( Intersector::polygon_area_2(n_vertices_result,  result_buffer) > 0 );
		
		intersection.resize(n_vertices_result, 2);
		std::copy(result_buffer, result_buffer + n_vertices_result * 2, &intersection.get_values()[0]);
		
		// plot_polygon(2, n_vertices_result, &intersection.get_values()[0], "isect/poly");
		return true;
	}
	
	void make_polygon_from_quad4(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		polygon.resize(e.n_nodes(), 2);
		
		for(int i = 0; i < e.n_nodes(); ++i) {
			for(int j = 0; j < 2; ++j) {
				polygon(i, j) = e.point(i)(j);
			}
		}
	}
	
	
	// void make_polygon_from_quad8(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	// {
	// 	polygon.resize(e.n_nodes()/2, 2);
		
	// 	for(int i = 0; i < e.n_nodes()/2; ++i) {
	// 		for(int j = 0; j < 2; ++j) {
	// 			polygon(i, j) = e.point(i)(j);
	// 			// std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
	// 		}
	// 	}
	// }

	void make_polygon_from_high_order_quad(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		polygon.resize(4, 2);
		
		for(int i = 0; i < 4; ++i) {
			for(int j = 0; j < 2; ++j) {
				polygon(i, j) = e.point(i)(j);
				// std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
			}
		}
	}
	
	void make_polygon_from_tri3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		polygon.resize(e.n_nodes(), 2);
		
		for(int i = 0; i < e.n_nodes(); ++i) {
			for(int j = 0; j < 2; ++j) {
				polygon(i, j) = e.point(i)(j);
			}
		}
	}
	
	template<int N>
	void discretize_segmented_curve(const double x[N], const double y[N],
									const int order,
									const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		using namespace libMesh;
		
		static const int Dim = 2;
		
		auto f = [&e](const double * x, double *fx) -> void {
			Point p;
			
			for(int d = 0; d < Dim; ++d) {
				p(d) = x[d];
			}
			
			Point fp = FE<Dim, libMesh::LAGRANGE>::map(&e, p);
			
			for(int d = 0; d < Dim; ++d) {
				fx[d] = fp(d);
			}
		};
		
		std::vector<double> all;
		std::vector<double> params_points, polyline;
		
		for(int i = 0; i < N; ++i) {
			const int ip1 = (i+1) % N;
			const double from[Dim] = { x[i],   y[i] };
			const double to  [Dim] = { x[ip1], y[ip1]};
			
			discretize_curve<Dim>(f, from, to, order, params_points, polyline, 1e-6);
			all.insert(all.end(), polyline.begin(), polyline.end()-2);
		}
		
		polygon.resize(all.size()/Dim, Dim);
		
		for(int i = 0; i < all.size(); i+= Dim) {
			polygon(i/Dim, 0) = all[i];
			polygon(i/Dim, 1) = all[i + 1];
		}
	}
	
	void make_polygon_from_curved_tri6(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		const double x[6] = { 0, 0.5, 1, 0.5, 0, 0.0 };
		const double y[6] = { 0, 0.0, 0, 0.5, 1, 0.5 };
		discretize_segmented_curve<6>(x, y, 2, e, polygon);
	}
	
	void make_polygon_from_curved_quad8(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		using namespace libMesh;
		const double x[8] = { -1,  0,  1, 1, 1, 0, -1, -1 };
		const double y[8] = { -1, -1, -1, 0, 1, 1,  1,  0 };
		discretize_segmented_curve<8>(x, y, 2, e, polygon);
	}
	
	void make_polygon_from_tri6(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		polygon.resize(e.n_nodes()/2, 2);
		for(int i = 0; i < e.n_nodes()/2; ++i) {
			for(int j = 0; j < 2; ++j) {
				polygon(i, j) = e.point(i)(j);
				// std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
				
			}
		}
	}
	
	void make_polygon(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		//FIXME use libMesh enum types
		switch(e.n_nodes()){
			case 2:
			{
				//works for lines too
				make_polygon_from_tri3(e, polygon);
				break;
			}
			case 3:
			{
				make_polygon_from_tri3(e, polygon);
				break;
			}
				
			case 4:
			{
				make_polygon_from_quad4(e, polygon);
				break;
			}
				
			case 6:
			{
				if(e.has_affine_map()) {
					make_polygon_from_tri6(e, polygon);
					// make_polygon_from_tri3(e, polygon);
				} else {
					make_polygon_from_curved_tri6(e, polygon);
				}
				break;
			}
				
			case 8:
			{
				if(e.has_affine_map()) {
					make_polygon_from_high_order_quad(e, polygon);
					// make_polygon_from_quad4(e, polygon);
				} else {
					make_polygon_from_curved_quad8(e, polygon);
				}
				break;
			}

			case 9:
			{
				// if(e.has_affine_map()) {
					make_polygon_from_high_order_quad(e, polygon);
					// make_polygon_from_quad4(e, polygon);
				// } else {
				// 	make_polygon_from_curved_quad8(e, polygon);
				// }
				break;
			}
				
			default:
			{
				assert(false);
				break;
			}
		}
	}
	
	void make_polygon_3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
	{
		auto n_nodes = e.n_nodes();

		if(e.has_affine_map()) {
			if(is_tri(e.type())) {
				n_nodes = 3;
			} else if(is_quad(e.type())) {
				n_nodes = 4;
			} else {
				assert(false && "handle special case");
			}
		} 

		polygon.resize(n_nodes, 3);
		
		for(int i = 0; i < n_nodes; ++i) {
			for(int j = 0; j < 3; ++j) {
				polygon(i, j) = e.point(i)(j);
			}
		}
	}


	static void make_polyhedron_from_generic_tet(const libMesh::Elem &e, Polyhedron &polyhedron)
	{		
		polyhedron.n_elements = 4;
		polyhedron.n_nodes 	  = 4;
		polyhedron.n_dims	  = 3;

		for(int i = 0; i < 4; ++i) {
			const int offset = i * 3;
			
			for(int j = 0; j < 3; ++j) {
				polyhedron.points[offset + j] = e.point(i)(j);
			}
		}
		
		polyhedron.el_ptr[0] = 0;
		polyhedron.el_ptr[1] = 3;
		polyhedron.el_ptr[2] = 6;
		polyhedron.el_ptr[3] = 9;
		polyhedron.el_ptr[4] = 12;
		
		//face 0
		polyhedron.el_index[0] = 0;
		polyhedron.el_index[1] = 1;
		polyhedron.el_index[2] = 3;
		
		//face 1
		polyhedron.el_index[3] = 1;
		polyhedron.el_index[4] = 2;
		polyhedron.el_index[5] = 3;
		
		//face 2
		polyhedron.el_index[6] = 0;
		polyhedron.el_index[7] = 3;
		polyhedron.el_index[8] = 2;
		
		//face 3
		polyhedron.el_index[9]  = 1;
		polyhedron.el_index[10] = 2;
		polyhedron.el_index[11] = 0;
	}
	
	void make_polyhedron_from_tet4(const libMesh::Elem &e, Polyhedron &polyhedron)
	{
		assert(e.dim() == 3);
		assert(e.n_nodes() == 4);

		make_polyhedron_from_generic_tet(e, polyhedron);
	}
	
	void make_polyhedron_from_tet10(const libMesh::Elem &e, Polyhedron &polyhedron)
	{
		assert(e.dim() == 3);
		assert(e.n_nodes() == 10);
		
		make_polyhedron_from_generic_tet(e, polyhedron);
	}

	static void make_polyhedron_from_generic_hex(const libMesh::Elem &e, Polyhedron &polyhedron)
	{
		polyhedron.n_elements = 6;
		polyhedron.n_nodes 	  = 8;
		polyhedron.n_dims	  = 3;
		
		
		for(int i = 0; i < 8; ++i) {
			const int offset = i * 3;
			
			for(int j = 0; j < 3; ++j) {
				polyhedron.points[offset + j] = e.point(i)(j);
			}
		}
		
		polyhedron.el_ptr[0] = 0;
		polyhedron.el_ptr[1] = 4;
		polyhedron.el_ptr[2] = 8;
		polyhedron.el_ptr[3] = 12;
		polyhedron.el_ptr[4] = 16;
		polyhedron.el_ptr[5] = 20;
		polyhedron.el_ptr[6] = 24;
		
		//face 0
		polyhedron.el_index[0] = 0;
		polyhedron.el_index[1] = 1;
		polyhedron.el_index[2] = 5;
		polyhedron.el_index[3] = 4;
		
		//face 1
		polyhedron.el_index[4] = 1;
		polyhedron.el_index[5] = 2;
		polyhedron.el_index[6] = 6;
		polyhedron.el_index[7] = 5;
		
		//face 2
		polyhedron.el_index[8]  = 3;
		polyhedron.el_index[9]  = 7;
		polyhedron.el_index[10] = 6;
		polyhedron.el_index[11] = 2;
		
		//face 3
		polyhedron.el_index[12] = 0;
		polyhedron.el_index[13] = 4;
		polyhedron.el_index[14] = 7;
		polyhedron.el_index[15] = 3;
		
		//face 4
		polyhedron.el_index[16] = 2;
		polyhedron.el_index[17] = 1;
		polyhedron.el_index[18] = 0;
		polyhedron.el_index[19] = 3;
		
		//face 5
		polyhedron.el_index[20] = 6;
		polyhedron.el_index[21] = 7;
		polyhedron.el_index[22] = 4;
		polyhedron.el_index[23] = 5;
	}
	
	
	//FIXME between -1, 1 ?
	void make_polyhedron_from_hex8(const libMesh::Elem &e, Polyhedron &polyhedron)
	{
		assert(e.dim() == 3);
		assert(e.n_nodes() == 8);
		
		make_polyhedron_from_generic_hex(e, polyhedron);
		
	}
	
	
	//FIXME between -1, 1 ?
	void make_polyhedron_from_hex27(const libMesh::Elem &e, Polyhedron &polyhedron)
	{
		assert(e.dim() == 3);
		assert(e.n_nodes() == 27);
		
		make_polyhedron_from_generic_hex(e, polyhedron);
	}
	
	void make_polyhedron(const libMesh::Elem &e, Polyhedron &polyhedron)
	{
		//FIXME use libMesh enum types
		switch(e.n_nodes()) {
			case 4:
			{
				make_polyhedron_from_tet4(e, polyhedron);
				break;
			}
				
			case 8:
			{
				make_polyhedron_from_hex8(e, polyhedron);
				break;
			}
				
			case 10:
			{
				make_polyhedron_from_tet10(e, polyhedron);
				break;
			}

			case 20:
			{
				make_polyhedron_from_generic_hex(e, polyhedron);
				break;
			}
				
			case 27:
			{
				make_polyhedron_from_hex27(e, polyhedron);
				break;
			}
				
			default:
			{
				assert(false);
				break;
			}
		}
	}
	
	bool intersect_3D(const libMesh::Elem &el1, const libMesh::Elem &el2, Polyhedron &intersection)
	{
		Polyhedron p1, p2;
		make_polyhedron(el1, p1);
		make_polyhedron(el2, p2);
		return Intersector::intersect_convex_polyhedra(p1, p2, &intersection);
	}
	
	bool intersect_3D(const Polyhedron &poly1, const Polyhedron &poly2, Polyhedron &intersection)
	{
		return Intersector::intersect_convex_polyhedra(poly1, poly2, &intersection);
	}
	
	bool project_2D(const libMesh::DenseMatrix<libMesh::Real> &poly1,
					const libMesh::DenseMatrix<libMesh::Real> &poly2,
					libMesh::DenseMatrix<libMesh::Real> &projection_1,
					libMesh::DenseMatrix<libMesh::Real> &projection_2)
	{
		using std::max;
		using std::min;
		
		typedef Intersector::Scalar Scalar;
		
		
		Scalar A   [2 * 2], b   [2];
		Scalar Ainv[2 * 2], binv[2];
		
//		Scalar normal_master[2];
//		Scalar normal_slave [2];
		
		Intersector::SurfaceMortarWorkspace w;
		
		const uint n_points_master = poly1.m();
		
		//////////////////////////////////////////////////////////////////////////////////////////
		//computing geometric surface projection
		
		//moving from global space to reference space
		Intersector::line_make_affine_transform_2(&poly2.get_values()[0], A, b);
		Intersector::make_inverse_affine_transform_2(A, b, Ainv, binv);
		
		Intersector::apply_affine_transform_2(Ainv, binv, n_points_master, &poly1.get_values()[0], w.ref_points_master);
		
		Scalar x_min, y_min;
		Scalar x_max, y_max;
		
		if(w.ref_points_master[0] < w.ref_points_master[2]) {
			x_min = w.ref_points_master[0];
			y_min = w.ref_points_master[1];
			
			x_max = w.ref_points_master[2];
			y_max = w.ref_points_master[3];
		} else {
			x_min = w.ref_points_master[2];
			y_min = w.ref_points_master[3];
			
			x_max = w.ref_points_master[0];
			y_max = w.ref_points_master[1];
		}
		
		//check if there is any intersection
		if(x_max <= 0) {
			return false;
		}
		
		if(x_min >= 1) {
			return false;
		}
		
		const Scalar x_min_isect = max(x_min, (Scalar)0);
		const Scalar x_max_isect = min(x_max, (Scalar)1);
		
		const Scalar dx = (x_max - x_min);
		const Scalar dy = (y_max - y_min);
		
		const Scalar y_min_isect = (x_min_isect - x_min) / dx * dy + y_min;
		const Scalar y_max_isect = (x_max_isect - x_min) / dx * dy + y_min;
		
		//store intersection lines
		w.isect_master[0] = x_min_isect;
		w.isect_master[1] = y_min_isect;
		w.isect_master[2] = x_max_isect;
		w.isect_master[3] = y_max_isect;
		
		w.isect_slave[0] = x_min_isect;
		w.isect_slave[1] = 0;
		w.isect_slave[2] = x_max_isect;
		w.isect_slave[3] = 0;
		
//		const Scalar inv_area_slave = 2.0/( Intersector::det_2(A) );
		
		//////////////////////////////////////////////////////////////////////////////////////////
		//create master fe object from intersection
		
		projection_1.resize(2, 2);
		projection_2.resize(2, 2);
		
		//move back to global coordinates
		Intersector::apply_affine_transform_2(A, b, 2, w.isect_master, &projection_1.get_values()[0]);
		
		//move back to global coordinates
		Intersector::apply_affine_transform_2(A, b, 2, w.isect_slave,  &projection_2.get_values()[0]);
		return true;
	}
	
	
	bool project_3D(const libMesh::DenseMatrix<libMesh::Real> &polygon_1,
					const libMesh::DenseMatrix<libMesh::Real> &polygon_2,
					libMesh::DenseMatrix<libMesh::Real> &projection_1,
					libMesh::DenseMatrix<libMesh::Real> &projection_2)
	{
		using namespace libMesh;
		
		typedef Intersector::Scalar Scalar;
		
		const int dim = 3;
		
		
		Scalar A   [3 * 3], b   [3];
		Scalar Ainv[3 * 3], binv[3];
		
//		Scalar normal_master[3];
//		Scalar normal_slave [3];
		
		Scalar isect_1[MAX_N_ISECT_POINTS * 3];
		Scalar isect_2[MAX_N_ISECT_POINTS * 3];
		
		libMesh::DenseMatrix<libMesh::Real> ref_polygon_1(polygon_1.m(), polygon_1.n());
		libMesh::DenseMatrix<libMesh::Real> ref_polygon_2(polygon_2.m(), polygon_2.n());
		libMesh::DenseMatrix<libMesh::Real> clipper_2(polygon_1.m(), 2);
		
		ref_polygon_1.resize(polygon_1.m(), polygon_1.n());
		clipper_2.resize(polygon_1.m(), dim - 1);
		ref_polygon_2.resize(polygon_2.m(), polygon_2.n());
		
		Intersector::triangle_make_affine_transform_3(&polygon_2.get_values()[0], A, b);
		Intersector::make_inverse_affine_transform_3(A, b, Ainv, binv);
		
		Intersector::apply_affine_transform_3(Ainv, binv,
										 polygon_1.m(),
										 &polygon_1.get_values()[0],
										 &ref_polygon_1.get_values()[0]);
		
		//FIXME could store reference instead of computing it each time
		Intersector::apply_affine_transform_3(Ainv, binv, polygon_2.m(),
										 &polygon_2.get_values()[0],
										 &ref_polygon_2.get_values()[0]);
		
		
		for(uint i = 0; i < polygon_1.m(); ++i) {
			clipper_2(i, 0) = ref_polygon_2(i, 0);
			clipper_2(i, 1) = ref_polygon_2(i, 1);
		}
		
		const int n_projection_points = Intersector::project_surface_poly_onto_ref_poly(
																				   dim,
																				   ref_polygon_1.m(), &ref_polygon_1.get_values()[0],
																				   polygon_2.m(),     &clipper_2.get_values()[0],
																				   isect_1, isect_2);
		
		
		// plot_polygon(clipper_2.n(),     clipper_2.m(),      &clipper_2.get_values()[0], "clipper");
		// plot_polygon(ref_polygon_1.n(), ref_polygon_1.m(),  &ref_polygon_1.get_values()[0], "clipped");
		// plot_polygon(3, n_projection_points,  isect_1, "clipper_isect");
		// plot_polygon(3, n_projection_points,  isect_2, "clipped_isect");
		
		if(!n_projection_points) {
			return false;
		}
		
		projection_1.resize(n_projection_points, 3);
		projection_2.resize(n_projection_points, 3);
		
		Intersector::apply_affine_transform_3(A, b, n_projection_points, isect_1, &projection_1.get_values()[0]);
		Intersector::apply_affine_transform_3(A, b, n_projection_points, isect_2, &projection_2.get_values()[0]);
		return true;
	}
	
	
	bool biorthgonal_weights(const int type, libMesh::Real &w_ii, libMesh::Real &w_ij)
	{
		using namespace libMesh;
		
		switch(type) {
			case EDGE2:
			{
				w_ii = 2.0;
				w_ij = -1.0;
				return true;
			}
				
			case TRI3:
			{
				w_ii = 3.0;
				w_ij = -1.0;
				return true;
			}
				
			case TET4:
			{
				w_ii = 4.0;
				w_ij = -1.0;
				return true;
			}
				// These do not work:
				// case QUAD4:
				// {
				// 	w_ii = 4.0;
				// 	w_ij = -1.0;
				// 	return true;
				// }
				
				// case HEX8:
				// {
				// 	w_ii = 8.0;
				// 	w_ij = -1.0;
				// 	return true;
				// }
				
			default:
			{
				w_ii = 1.0;
				w_ij = 0.0;
				
				static bool error_msg_printed = false;
				
				if(!error_msg_printed) {
					std::cerr << "[Error] biorthgonal weights not supported for element type: " << type << std::endl;
					error_msg_printed = true;
				}
				
				assert(false && "TODO: add the weights for the missing element");
				return false;
			}
		}
	}
	
	template<class FE>
	void mortar_assemble_weights_aux(const FE &fe, libMesh::DenseMatrix<libMesh::Real> &weights)
	{
		libMesh::DenseMatrix<libMesh::Real> elmat;
		elmat.resize(fe.get_phi().size(), fe.get_phi().size());
		elmat.zero();
		
		
		weights.resize(elmat.m(), elmat.n());
		weights.zero();
		
		const auto &test = fe.get_phi();
		const auto &JxW   = fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_qp    = test[0].size();
		
		for(uint qp = 0; qp < n_qp; ++qp) {
			for(uint i = 0; i < n_test; ++i) {
				for(uint j = 0; j < n_test; ++j) {
					elmat(i, j) += contract(test[i][qp], test[j][qp]) * JxW[qp];
				}
			}
		}
		
		libMesh::DenseVector<libMesh::Real> sum_elmat(n_test);
		sum_elmat.zero();
		libMesh::DenseVector<libMesh::Real> rhs(n_test);
		rhs.zero();
		
		libMesh::DenseVector<libMesh::Real> sol(n_test);
		sol.zero();
		
		for(uint i = 0; i < n_test; ++i) {
			for(uint j = 0; j < n_test; ++j) {
				sum_elmat(i) += elmat(i, j);
			}
			
			if(std::abs(sum_elmat(i)) < 1e-16) {
				sum_elmat(i) = 0;
				//set identity row where not defined
				for(uint j = 0; j < n_test; ++j) {
					elmat(i, j) = (i == j);
				}
			}
		}
		
		// std::cout << "-----------------------\n";
		// std::cout << "-----------------------\n";
		
		// elmat.print(std::cout);
		
		// std::cout << "-----------------------\n";
		
		for(uint i = 0; i < n_test; ++i) {
			if(sum_elmat(i) == 0) {
				continue;
			}
			
			rhs(i) = sum_elmat(i);
			
			elmat.cholesky_solve(rhs, sol);
			
			for(uint j = 0; j < n_test; ++j) {
				weights(i, j) = sol(j);
			}
			
			rhs(i) = 0;
		}
		
		//normalization for consistently scaled coefficients
		for(uint i = 0; i < n_test; ++i) {
			if(sum_elmat(i) == 0) {
				continue;
			}
			
			libMesh::Real t = 0;
			for(uint j = 0; j < n_test; ++j) {
				t += weights(i, j);
			}
			
			for(uint j = 0; j < n_test; ++j) {
				weights(i, j) *= 1./t;
			}
		}
		
		// weights.print(std::cout);
	}
	
	void mortar_assemble_weights(const libMesh::FEVectorBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights)
	{
		mortar_assemble_weights_aux(fe, weights);
	}
	
	void mortar_assemble_weights(const libMesh::FEBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights)
	{
		mortar_assemble_weights_aux(fe, weights);
	}
	
	
	template<class FE>
	void mortar_assemble_weighted_aux(
									  const FE &trial_fe,
									  const FE &test_fe,
									  const libMesh::DenseMatrix<libMesh::Real> &weights,
									  libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		if(elmat.m() != test_fe.get_phi().size() ||
		   elmat.n() != trial_fe.get_phi().size()) {
			
			elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
			elmat.zero();
		}
		
		const auto &trial = trial_fe.get_phi();
		const auto &test  = test_fe.get_phi();
		const auto &JxW   = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_trial = trial.size();
		const uint n_qp    = test[0].size();
		
		for(uint i = 0; i < n_test; ++i) {
			for(uint qp = 0; qp < n_qp; ++qp) {
				auto w_test = test[i][qp] * 0.;
				
				for(uint k = 0; k < n_test; ++k) {
					w_test += test[k][qp] * weights(i, k);
				}
				
				for(uint j = 0; j < n_trial; ++j) {
					// assert(  JxW[qp] >= 0.0 );
					elmat(i, j) += contract(w_test, trial[j][qp]) * JxW[qp];
				}
			}
		}
	}
	
	void mortar_assemble_weighted_biorth(
										 const libMesh::FEBase &trial_fe,
										 const libMesh::FEBase &test_fe,
										 const libMesh::DenseMatrix<libMesh::Real> &weights,
										 libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		mortar_assemble_weighted_aux(trial_fe, test_fe, weights, elmat);
	}
	
	
	
	
	void mortar_assemble_weighted_biorth(
										 const libMesh::FEVectorBase &trial_fe,
										 const libMesh::FEVectorBase &test_fe,
										 const libMesh::DenseMatrix<libMesh::Real> &weights,
										 libMesh::DenseMatrix<libMesh::Real> &elmat)
	{
		mortar_assemble_weighted_aux(trial_fe, test_fe, weights, elmat);
	}
	
	void mortar_normal_and_gap_assemble_weighted_biorth(
														const libMesh::FEVectorBase &test_fe,
														const int dim,
														const libMesh::Point &surf_normal,
														const libMesh::Point &plane_normal,
														const libMesh::Real &plane_offset,
														const libMesh::DenseMatrix<libMesh::Real> &weights,
														libMesh::DenseMatrix<libMesh::Real> &normals,
														libMesh::DenseVector<libMesh::Real> &gap)
	{
		
		using namespace libMesh;
		
		
		
		if(normals.m() != test_fe.get_phi().size()/dim || dim != normals.n()) {
			normals.resize(test_fe.get_phi().size()/dim, dim);
			gap.resize(test_fe.get_phi().size());
		}
		
		normals.zero();
		gap.zero();
		
		const auto &test   = test_fe.get_phi();
		const auto &point  = test_fe.get_xyz();
		const auto &JxW    = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_qp    = test[0].size();
		
		DenseVector<Real> p(dim);
		DenseVector<Real> v(dim);
		
		DenseVector<Real> s_n(dim);
		DenseVector<Real> p_n(dim);
		
		for(uint i = 0; i < dim; ++i) {
			p_n(i) =  plane_normal(i);
			s_n(i) =  surf_normal(i);
		}
		
		for(uint qp = 0; qp < n_qp; ++qp) {
			
			p(0) = point[qp](0);
			p(1) = point[qp](1);
			
			if(dim > 2) {
				p(2) = point[qp](2);
			}
			
			Real isect = 0;
			Intersector::intersect_ray_with_plane(dim, 1, &p.get_values()[0], &s_n.get_values()[0], &p_n.get_values()[0], plane_offset, &isect);
			
			v = s_n;
			v *= isect;
			// quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);
			
			for(uint i = 0; i < n_test; ++i) {
				auto biorth_test =  weights(i, 0) * test[0][qp];
				
				for(uint k = 0; k < test.size(); ++k) {
					biorth_test +=  weights(i, k) * test[k][qp];
				}
				
				gap(i) += biorth_test(0) * isect * JxW[qp];
				
				for(uint d = 0; d < dim; ++d) {
					normals.get_values()[i] += biorth_test(d) * surf_normal(d) * JxW[qp];
				}
			}
		}
		
		// gap.print(std::cout);
		
	}
	
	template<typename T>
	void convert_point_to_vector(const int dim, const libMesh::Point &point, std::vector<T> &point_vec)
	{
		point_vec.resize(dim);
		for(int i = 0; i < dim; ++i) {
			point_vec[i] = point(i);
		}
	}
	
	void mortar_normal_and_gap_assemble_weighted_biorth(
														const libMesh::FEBase &test_fe,
														const int dim,
														const libMesh::Point &surf_normal,
														const libMesh::Point &plane_normal,
														const libMesh::Real &plane_offset,
														const libMesh::DenseMatrix<libMesh::Real> &weights,
														libMesh::DenseMatrix<libMesh::Real> &normals,
														libMesh::DenseVector<libMesh::Real> &gap,
														const bool visdbg
														)
	{
		using namespace libMesh;
		
		
		if(normals.m() != test_fe.get_phi().size() || dim != normals.n()) {
			normals.resize(test_fe.get_phi().size(), dim);
			normals.zero();
			gap.resize(test_fe.get_phi().size());
			gap.zero();
		}
		
		const auto &test   = test_fe.get_phi();
		// const auto &grad   = test_fe.get_dphi();
		const auto &point  = test_fe.get_xyz();
		const auto &JxW    = test_fe.get_JxW();
		
		const uint n_test  = test.size();
		const uint n_qp    = test[0].size();
		
		
		std::vector<Real> surf_normal_v, plane_normal_v;
		convert_point_to_vector(dim, surf_normal, surf_normal_v);
		convert_point_to_vector(dim, plane_normal, plane_normal_v);
		
		DenseVector<Real> p(dim);
		// DenseMatrix<Real> v(dim); //visdbg
		
		for(uint i = 0; i < n_test; ++i) {
			for(uint qp = 0; qp < n_qp; ++qp) {
				
				p(0) = point[qp](0);
				p(1) = point[qp](1);
				
				if(dim > 2) {
					p(2) = point[qp](2);
				}
				
				Real isect = 0;
				Intersector::intersect_ray_with_plane(dim, 1, &p.get_values()[0], &surf_normal_v[0], &plane_normal_v[0], plane_offset, &isect);
				// assert(isect > 0);
				// if(visdbg) {
				// 	v.get_values() = surf_normal_v; //visdbg
				// 	v *= isect; //visdbg
				// 	quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]); //visdbg
				// }
				
				auto biorth_test = weights(i, 0) * test[0][qp];
				
				for(uint k = 1; k < n_test; ++k) {
					biorth_test += weights(i, k)  * test[k][qp];
				}
				
				
				gap(i) += biorth_test * isect * JxW[qp];
				
				for(uint d = 0; d < dim; ++d) {
					normals(i, d) += biorth_test * surf_normal(d) * JxW[qp];
				}
			}
		}
	}
}
