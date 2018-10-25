#include "utopia_ApproxL2LocalAssembler.hpp"

namespace utopia {

	using Vector2 = Intersector::Vector2;
	inline static bool inside_half_plane(const Vector2 &e1, const Vector2 &e2, const Vector2 &point, const double tol)
	{
		const Vector2 u = e1 - e2;
		const Vector2 v = point - e2;

		const double dist = (u.x * v.y) - (v.x * u.y);
		return dist <= tol;
	}

	std::shared_ptr<Transform> ApproxL2LocalAssembler::get_trafo(const Elem &elem) const
	{
		std::shared_ptr<Transform> elem_trafo;
		//FIXME
		// if(elem.has_affine_map()) {
		// 	if(dim == 2) {
		// 		if(elem.dim() == 1) {
		// 			elem_trafo = std::make_shared<Transform1>(elem);
		// 		} else {
		// 			elem_trafo = std::make_shared<AffineTransform2>(elem);
		// 		}
		// 	} else {
		// 		assert(dim == 3);

		// 		if(elem.dim() == 2) {
		// 			elem_trafo = std::make_shared<Transform2>(elem);
		// 		} else {
		// 			elem_trafo = std::make_shared<AffineTransform3>(elem);
		// 		}
		// 	}
		// } else {
			if(elem.dim() == 1) {
				elem_trafo = std::make_shared<Transform1>(elem);
			} else if(dim == 2) {
				assert(elem.dim() == 2);
				elem_trafo = std::make_shared<Transform2>(elem);
			} else {
				assert(dim == 3);

				if(elem.dim() == 2) {
					elem_trafo = std::make_shared<Transform2>(elem);
				} else {
					elem_trafo = std::make_shared<Transform3>(elem);
				}
			}
		// }

		return elem_trafo;
	}

	bool ApproxL2LocalAssembler::assemble(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		Matrix &mat
		)
	{
		if(!init_q(trial, trial_type, test, test_type)) {
			return false;
		}

		auto trial_fe = libMesh::FEBase::build(trial.dim(), trial_type);
		trial_fe->attach_quadrature_rule(q_trial.get());
		trial_fe->get_phi();
		trial_fe->reinit(&trial);

		auto test_fe = libMesh::FEBase::build(test.dim(), test_type);
		test_fe->attach_quadrature_rule(q_test.get());
		test_fe->get_JxW();
		test_fe->get_phi();
		test_fe->reinit(&test);

		assemble(*trial_fe, *test_fe, mat);

		assert((test_type != libMesh::FIRST || trial_type != libMesh::FIRST || check_valid(mat)));
		return true;
	}

	bool ApproxL2LocalAssembler::assemble(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		std::vector<Matrix> &mat
		) 
	{
		assert(n_forms() == 2);

		mat.resize(n_forms());

		if(!init_q(trial, trial_type, test, test_type)) {
			return false;
		}

		auto trial_fe = libMesh::FEBase::build(trial.dim(), trial_type);
		trial_fe->attach_quadrature_rule(q_trial.get());
		trial_fe->get_phi();
		trial_fe->reinit(&trial);

		auto test_fe = libMesh::FEBase::build(test.dim(), test_type);
		test_fe->attach_quadrature_rule(q_test.get());
		test_fe->get_JxW();
		test_fe->get_phi();
		test_fe->reinit(&test);

		assemble(*trial_fe, *test_fe, mat[0]);
		assemble(*test_fe,  *test_fe, mat[1]);
		return true;
	}

	void ApproxL2LocalAssembler::assemble(libMesh::FEBase &trial, libMesh::FEBase &test, Matrix &mat) const
	{
		const auto &trial_shape_fun = trial.get_phi();
		const auto &test_shape_fun  = test.get_phi();
		auto &JxW = test.get_JxW();

		mat.resize(test_shape_fun.size(), trial_shape_fun.size());

		std::size_t n_qps = JxW.size();
		for(std::size_t i = 0; i < test_shape_fun.size(); ++i) {
			for(std::size_t j = 0; j < trial_shape_fun.size(); ++j) {
				for(std::size_t k = 0; k < n_qps; ++k) {
					auto tf = test_shape_fun.at(i).at(k);
					mat(i, j) += trial_shape_fun.at(j).at(k) * tf * JxW[k];
				}
			}
		}
	}

	bool ApproxL2LocalAssembler::check_valid(const Matrix &mat) const
	{
		for(int i = 0; i < mat.m(); ++i) {
			double row_sum = 0.;
			for(int j = 0; j < mat.n(); ++j) {
				row_sum += mat(i, j);
				assert(mat(i, j) < 1.0001);
			}

			assert(row_sum < 1.0001);
			if(row_sum > 1.001) return false;
		}

		return true;
	}

	bool ApproxL2LocalAssembler::init_q(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type)
	{
		std::shared_ptr<Transform> trial_trafo = get_trafo(trial);
		std::shared_ptr<Transform> test_trafo  = get_trafo(test);

		if(!q_trial) {
			q_trial = std::make_shared<QMortar>(dim);
			q_test  = std::make_shared<QMortar>(dim);
		}

		int order = quadrature_order;
		if(order < 0) {
			order = std::max(
				order_for_l2_integral(dim, trial, trial_type.order, test, test_type.order),
				order_for_l2_integral(dim, test,  test_type.order,  test, test_type.order)
			);
		}

		libMesh::QGauss q(test.dim(), libMesh::Order(order));
		q.init(test.type());

		std::size_t n_quad_points = q.get_points().size();

		QMortar q_spatial(dim);
		q_spatial.resize(n_quad_points);

		for(std::size_t i = 0; i < n_quad_points; ++i) {
			test_trafo->apply(q.get_points()[i], q_spatial.get_points()[i]);
		}

		std::vector<int> index;
		contained_points(trial, q_spatial, index);

		if(index.empty()) return false;

		std::size_t n_isect_qp = index.size();

		q_trial->resize(n_isect_qp);
		q_test->resize(n_isect_qp);

		std::fill(q_trial->get_weights().begin(), q_trial->get_weights().end(), 0.);
		std::fill(q_test->get_weights().begin(),  q_test->get_weights().end(), 0.);

		for(std::size_t qp = 0; qp < n_isect_qp; ++qp) {
			auto ind = index[qp];
			auto p = q_spatial.get_points()[ind];

			trial_trafo->transform_to_reference(p, q_trial->get_points()[qp]);
			q_test->get_points()[qp]  = q.get_points()[ind];
			q_test->get_weights()[qp] = q.get_weights()[ind];
		}

		return true;
	}

	void ApproxL2LocalAssembler::contained_points(const Elem &trial, const libMesh::QBase &q, std::vector<int> &index)
	{
		int n_potential_nodes = q.get_points().size();

		// if(nested_meshes) {
		// 	//check if there is an intersection then...


		// 	test_dofs.resize(n_potential_nodes);

		// 	for(int i = 0; i < n_potential_nodes; ++i) {
		// 		test_dofs[i] = i;
		// 	}

		// 	return;
		// }

		index.reserve(n_potential_nodes);

		if(dim == 2) {
			contained_points_2(trial, q, index);
			return;
		} else if(dim == 3) {
			assert(dim == 3);
			contained_points_3(trial, q, index);
			return;
		}

		for(int i = 0; i < n_potential_nodes; ++i) {
			auto const & q_node = q.get_points()[i];
			if(trial.contains_point(q_node, tol)) {
				index.push_back(i);
			}
		}
	}

	void ApproxL2LocalAssembler::contained_points_2(const Elem &trial, const libMesh::QBase &q, std::vector<int> &index)
	{
		make_polygon(trial, trial_pts);
		// make_polygon(test,  test_pts);

		const auto &trial_poly = trial_pts.get_values();
		const int trial_n_nodes = trial_poly.size() / 2;
		const int n_qp  = q.get_points().size();

		index.clear();

		Vector2 e1, e2, p, s, e;

		std::vector<bool> is_inside(n_qp, true);

		for(int i = 0; i < trial_n_nodes; ++i) {
			const int i2x = i * 2;
			const int i2y = i2x + 1;

			const int i2p1x = 2 * (((i + 1) == trial_n_nodes)? 0 : (i + 1));
			const int i2p1y = i2p1x + 1;

			e1 = Vector2(trial_poly[i2x],   trial_poly[i2y]);
			e2 = Vector2(trial_poly[i2p1x], trial_poly[i2p1y]);

			for(int j = 0; j < n_qp; ++j) {
				const auto &q_pt = q.get_points()[j];
				p = Vector2(q_pt(0), q_pt(1));
				is_inside[j] = is_inside[j] && inside_half_plane(e1, e2, p, tol);
			}
		}

		for(int i = 0; i < n_qp; ++i) {
			if(is_inside[i]) {
				index.push_back(i);
			}
		}
	}

	void ApproxL2LocalAssembler::contained_points_3(const Elem &trial,  const libMesh::QBase &q, std::vector<int> &index) const
	{
		Polyhedron poly;
		make_polyhedron(trial, poly);

		auto n_half_spaces = poly.n_elements;

		std::vector<double> plane_normals(n_half_spaces * 3, 0.);
		std::vector<double> plane_dists_from_origin(n_half_spaces, 0.);
		Intersector::make_h_polyhedron_from_polyhedron(poly, &plane_normals[0], &plane_dists_from_origin[0]);

		double p[3];

		index.clear();

		int n_qp = q.get_points().size();
		for(int k = 0; k < n_qp; ++k) {
			const auto & node = q.get_points()[k];
			p[0] = node(0);
			p[1] = node(1);
			p[2] = node(2);

			bool inside = true;
			for(int i = 0; i < n_half_spaces; ++i) {
				const int i3 = i * 3;
				auto d = Intersector::point_plane_distance(3, &plane_normals[i3], plane_dists_from_origin[i], p);

				if(d >= tol) {
					inside = false;
					break;
				}
			}

			if(inside) {
				index.push_back(k);
			}
		}
	}
}
