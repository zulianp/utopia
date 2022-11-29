// #include "MortarAssemble.hpp"
// #include <libmesh/fe.h>
// #include "utopia_Polygon.hpp"
// #include "utopia_libmesh_Deprecated.hpp"
// #include "utopia_libmesh_Transform.hpp"
// #include "utopia_libmesh_Utils.hpp"
// #include "utopia_triangulate.hpp"

// #include "Box.hpp"

// #include <assert.h>
// #include <algorithm>
// #include <memory>
// #include <numeric>

// #include "utopia_intersector.hpp"

// namespace utopia {

//     double ref_volume(int type) {
//         if (is_hex(type)) {
//             return 8.;
//         } else if (is_tet(type)) {
//             return 1. / 6.;
//         } else if (is_quad(type)) {
//             return 4.;
//         } else if (is_edge(type)) {
//             return 2.;
//         } else if (is_tri(type)) {
//             return 0.5;
//         } else if (is_prism(type)) {
//             return 1.;
//         } else if (is_pyramid(type)) {
//             return 1. / 0.75;
//         } else {
//             assert(false);
//             return 1.;
//         }
//     }

//     double ref_area_of_surf(int type) {
//         if (is_tri(type)) {
//             return 2.;
//         } else if (is_quad(type)) {
//             return 2.;
//         } else if (is_hex(type)) {
//             return 1.;
//         } else if (is_tet(type)) {
//             return 0.5;
//         } else if (is_edge(type)) {
//             // Is this correct?
//             assert(false);
//             return 0.5;
//         } else {
//             assert(false && "add special case");
//             return 1.;
//         }
//     }

//     void transform_to_reference(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir) {
//         const int dim = global_ir.get_dim();
//         libMesh::Point p;
//         ref_ir.resize(global_ir.n_points());

//         // libMesh::DenseMatrix<libMesh::Real> A_inv;

//         const double factor = ref_volume(type);

//         double sum_of_weights = 0.;
//         for (int i = 0; i < global_ir.n_points(); ++i) {
//             p(0) = global_ir.qp(i)(0);
//             p(1) = global_ir.qp(i)(1);
//             p(2) = global_ir.qp(i)(2);

//             trans.transform_to_reference(p, ref_ir.get_points()[i]);

//             ref_ir.get_weights()[i] = global_ir.w(i) * factor;

//             // for debugging
//             sum_of_weights += ref_ir.get_weights()[i];
//         }

//         // assert(sum_of_weights > 0.);
//         // assert(sum_of_weights <= factor + 5e-4);
//     }

//     void transform_to_reference_surf(const Transform &trans,
//                                      const int type,
//                                      const QMortar &global_ir,
//                                      QMortar &ref_ir) {
//         assert(is_valid_elem_type(type));

//         const int dim = global_ir.get_dim();
//         libMesh::Point p;
//         ref_ir.resize(global_ir.n_points());

//         libMesh::DenseMatrix<libMesh::Real> A_inv;

//         const double factor = ref_area_of_surf(type);

//         for (int i = 0; i < global_ir.n_points(); ++i) {
//             p(0) = global_ir.qp(i)(0);
//             p(1) = global_ir.qp(i)(1);
//             p(2) = global_ir.qp(i)(2);

//             trans.transform_to_reference(p, ref_ir.get_points()[i]);
//             ref_ir.get_weights()[i] = global_ir.w(i) * factor;
//         }
//     }

//     template <class Left, class Right>
//     libMesh::Real contract(const Left &left, const Right &right) {
//         return left.contract(right);
//     }

//     libMesh::Real contract(const libMesh::Real &left, const libMesh::Real &right) { return left * right; }

//     template <class FE>
//     void mortar_assemble_aux(const FE &trial_fe, const FE &test_fe, libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
//             elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
//             elmat.zero();
//         }

//         const auto &trial = trial_fe.get_phi();
//         const auto &test = test_fe.get_phi();
//         const auto &JxW = test_fe.get_JxW();

//         const uint n_test = test.size();
//         const uint n_trial = trial.size();
//         const uint n_qp = test[0].size();

//         assert(test[0].size() == trial[0].size());

//         for (uint i = 0; i < n_test; ++i) {
//             for (uint j = 0; j < n_trial; ++j) {
//                 for (uint qp = 0; qp < n_qp; ++qp) {
//                     elmat(i, j) += contract(test[i][qp], trial[j][qp]) * JxW[qp];
//                 }
//             }
//         }
//     }

//     libMesh::Real len(const libMesh::Real val) { return std::abs(val); }

//     template <class Vec>
//     libMesh::Real len(const Vec &val) {
//         return val.size();
//     }

//     static inline bool is_vec(const libMesh::Real val) { return false; }

//     template <class Vec>
//     bool is_vec(const Vec &val) {
//         return true;
//     }

//     void make_tp(const int i, libMesh::Real &val) {}

//     template <class Vec>
//     void make_tp(const int i, Vec &val) {
//         libMesh::Real s = val(i);
//         val.zero();
//         val(i) = s;
//     }

//     template <class FE>
//     void mortar_assemble_biorth_aux(const FE &trial_fe,
//                                     const FE &test_fe,
//                                     const libMesh::Real &wii,
//                                     const libMesh::Real &wij,
//                                     libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
//             elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
//             elmat.zero();
//         }

//         const auto &trial = trial_fe.get_phi();
//         const auto &test = test_fe.get_phi();
//         const auto &JxW = test_fe.get_JxW();

//         const uint n_test = test.size();
//         const uint n_trial = trial.size();
//         const uint n_qp = test[0].size();

//         for (uint qp = 0; qp < n_qp; ++qp) {
//             for (uint i = 0; i < n_test; ++i) {
//                 auto biorth_test = ((0 == i) ? wii : wij) * test[0][qp];

//                 for (uint k = 1; k < n_test; ++k) {
//                     biorth_test += ((k == i) ? wii : wij) * test[k][qp];
//                 }

//                 for (uint j = 0; j < n_trial; ++j) {
//                     elmat(i, j) += contract(biorth_test, trial[j][qp]) * JxW[qp];
//                 }
//             }
//         }
//     }

//     void mortar_assemble(const libMesh::FEBase &trial_fe,
//                          const libMesh::FEBase &test_fe,
//                          libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         mortar_assemble_aux(trial_fe, test_fe, elmat);
//     }

//     void mortar_assemble(const libMesh::FEVectorBase &trial_fe,
//                          const libMesh::FEVectorBase &test_fe,
//                          libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         mortar_assemble_aux(trial_fe, test_fe, elmat);
//     }

//     void mortar_assemble_biorth(const libMesh::FEBase &trial_fe,
//                                 const libMesh::FEBase &test_fe,
//                                 const int type,
//                                 libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         libMesh::Real w_ii, w_ij;
//         biorthgonal_weights(type, w_ii, w_ij);
//         mortar_assemble_biorth_aux(trial_fe, test_fe, w_ii, w_ij, elmat);
//     }

//     void mortar_assemble_biorth(const libMesh::FEVectorBase &trial_fe,
//                                 const libMesh::FEVectorBase &test_fe,
//                                 const int type,
//                                 libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         libMesh::Real w_ii, w_ij;
//         biorthgonal_weights(type, w_ii, w_ij);
//         mortar_assemble_biorth_aux(trial_fe, test_fe, w_ii, w_ij, elmat);
//     }

//     template <class FE>
//     void mortar_assemble_biorth_aux(const int dim,
//                                     const FE &trial_fe,
//                                     const FE &test_fe,
//                                     const libMesh::Real &wii,
//                                     const libMesh::Real &wij,
//                                     const libMesh::DenseVector<libMesh::Real> &indicator,
//                                     libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
//             elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
//             elmat.zero();
//         }

//         const auto &trial = trial_fe.get_phi();
//         const auto &test = test_fe.get_phi();
//         const auto &JxW = test_fe.get_JxW();

//         const uint n_test = test.size();
//         const uint n_trial = trial.size();
//         const uint n_qp = test[0].size();

//         //		bool v = is_vec(test[0][0]);

//         for (uint qp = 0; qp < n_qp; ++qp) {
//             for (uint i = 0; i < n_test; ++i) {
//                 auto biorth_test = ((0 == i) ? wii : wij) * test[0][qp];

//                 for (uint k = 1; k < n_test; ++k) {
//                     biorth_test += ((k == i) ? wii : wij) * test[k][qp];
//                 }

//                 for (int k = 0; k < dim; ++k) {
//                     if (i % dim == k) {
//                         make_tp(k, biorth_test);
//                         break;
//                     }
//                 }

//                 // if(indicator(i) > 0)
//                 // 	std::cout <<  biorth_test << std::endl;

//                 for (uint j = 0; j < n_trial; ++j) {
//                     elmat(i, j) += indicator(i) * contract(biorth_test, trial[j][qp]) * JxW[qp];
//                 }
//             }
//         }
//     }

//     void mortar_assemble_biorth(const int dim,
//                                 const libMesh::FEBase &trial_fe,
//                                 const libMesh::FEBase &test_fe,
//                                 const int type,
//                                 const libMesh::DenseVector<libMesh::Real> &indicator,
//                                 libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         libMesh::Real w_ii, w_ij;
//         biorthgonal_weights(type, w_ii, w_ij);
//         mortar_assemble_biorth_aux(dim, trial_fe, test_fe, w_ii, w_ij, indicator, elmat);
//     }

//     void mortar_assemble_biorth(const int dim,
//                                 const libMesh::FEVectorBase &trial_fe,
//                                 const libMesh::FEVectorBase &test_fe,
//                                 const int type,
//                                 const libMesh::DenseVector<libMesh::Real> &indicator,
//                                 libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         libMesh::Real w_ii, w_ij;
//         biorthgonal_weights(type, w_ii, w_ij);
//         mortar_assemble_biorth_aux(dim, trial_fe, test_fe, w_ii, w_ij, indicator, elmat);
//     }

//     void mortar_normal_and_gap_assemble_biorth(const int type,
//                                                const uint dim,
//                                                const libMesh::FEBase &test_fe,
//                                                const libMesh::Point &surf_normal,
//                                                const libMesh::Point &plane_normal,
//                                                const libMesh::Real &plane_offset,
//                                                const libMesh::DenseVector<libMesh::Real> &indicator,
//                                                libMesh::DenseMatrix<libMesh::Real> &normals,
//                                                libMesh::DenseVector<libMesh::Real> &gap) {
//         using namespace libMesh;
//         DenseVector<Real> surf_normal_v(dim), plane_normal_v(dim);

//         for (uint i = 0; i < dim; ++i) {
//             surf_normal_v(i) = surf_normal(i);
//             plane_normal_v(i) = plane_normal(i);
//         }

//         mortar_normal_and_gap_assemble_biorth(
//             type, test_fe, surf_normal_v, plane_normal_v, plane_offset, indicator, normals, gap);
//     }

//     void make_polygon_from_quad4(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         polygon.resize(e.n_nodes(), 2);

//         for (int i = 0; i < e.n_nodes(); ++i) {
//             for (int j = 0; j < 2; ++j) {
//                 polygon(i, j) = e.point(i)(j);
//             }
//         }
//     }

//     // void make_polygon_from_quad8(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon)
//     // {
//     // 	polygon.resize(e.n_nodes()/2, 2);

//     // 	for(int i = 0; i < e.n_nodes()/2; ++i) {
//     // 		for(int j = 0; j < 2; ++j) {
//     // 			polygon(i, j) = e.point(i)(j);
//     // 			// std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
//     // 		}
//     // 	}
//     // }

//     void make_polygon_from_high_order_quad(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         polygon.resize(4, 2);

//         for (int i = 0; i < 4; ++i) {
//             for (int j = 0; j < 2; ++j) {
//                 polygon(i, j) = e.point(i)(j);
//                 // std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
//             }
//         }
//     }

//     void make_polygon_from_tri3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         polygon.resize(e.n_nodes(), 2);

//         for (int i = 0; i < e.n_nodes(); ++i) {
//             for (int j = 0; j < 2; ++j) {
//                 polygon(i, j) = e.point(i)(j);
//             }
//         }
//     }

//     template <int N>
//     void discretize_segmented_curve(const double x[N],
//                                     const double y[N],
//                                     const int order,
//                                     const libMesh::Elem &e,
//                                     libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         using namespace libMesh;

//         static const int Dim = 2;

//         auto f = [&e](const double *x, double *fx) -> void {
//             Point p;

//             for (int d = 0; d < Dim; ++d) {
//                 p(d) = x[d];
//             }

//             Point fp = FE<Dim, libMesh::LAGRANGE>::map(&e, p);

//             for (int d = 0; d < Dim; ++d) {
//                 fx[d] = fp(d);
//             }
//         };

//         std::vector<double> all;
//         std::vector<double> params_points, polyline;

//         for (int i = 0; i < N; ++i) {
//             const int ip1 = (i + 1) % N;
//             const double from[Dim] = {x[i], y[i]};
//             const double to[Dim] = {x[ip1], y[ip1]};

//             discretize_curve<Dim>(f, from, to, order, params_points, polyline, 1e-6);
//             all.insert(all.end(), polyline.begin(), polyline.end() - 2);
//         }

//         polygon.resize(all.size() / Dim, Dim);

//         for (int i = 0; i < all.size(); i += Dim) {
//             polygon(i / Dim, 0) = all[i];
//             polygon(i / Dim, 1) = all[i + 1];
//         }
//     }

//     void make_polygon_from_curved_tri6(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         const double x[6] = {0, 0.5, 1, 0.5, 0, 0.0};
//         const double y[6] = {0, 0.0, 0, 0.5, 1, 0.5};
//         discretize_segmented_curve<6>(x, y, 2, e, polygon);
//     }

//     void make_polygon_from_curved_quad8(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         using namespace libMesh;
//         const double x[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
//         const double y[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
//         discretize_segmented_curve<8>(x, y, 2, e, polygon);
//     }

//     void make_polygon_from_tri6(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         polygon.resize(e.n_nodes() / 2, 2);
//         for (int i = 0; i < e.n_nodes() / 2; ++i) {
//             for (int j = 0; j < 2; ++j) {
//                 polygon(i, j) = e.point(i)(j);
//                 // std::cout<<" polygon("<<i<<","<<j<<") = "<< polygon(i, j) <<std::endl;
//             }
//         }
//     }

//     void make_polyline(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polyline) {
//         polyline.resize(2, 2);

//         polyline(0, 0) = e.node_ref(0)(0);
//         polyline(0, 1) = e.node_ref(0)(1);

//         polyline(1, 0) = e.node_ref(1)(0);
//         polyline(1, 1) = e.node_ref(1)(1);
//     }

//     void make_polygon(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         // FIXME use libMesh enum types
//         switch (e.n_nodes()) {
//             case 2: {
//                 // works for lines too
//                 make_polygon_from_tri3(e, polygon);
//                 break;
//             }
//             case 3: {
//                 make_polygon_from_tri3(e, polygon);
//                 break;
//             }

//             case 4: {
//                 make_polygon_from_quad4(e, polygon);
//                 break;
//             }

//             case 6: {
//                 if (e.has_affine_map()) {
//                     make_polygon_from_tri6(e, polygon);
//                     // make_polygon_from_tri3(e, polygon);
//                 } else {
//                     make_polygon_from_curved_tri6(e, polygon);
//                 }
//                 break;
//             }

//             case 8: {
//                 if (e.has_affine_map()) {
//                     make_polygon_from_high_order_quad(e, polygon);
//                     // make_polygon_from_quad4(e, polygon);
//                 } else {
//                     make_polygon_from_curved_quad8(e, polygon);
//                 }
//                 break;
//             }

//             case 9: {
//                 // if(e.has_affine_map()) {
//                 make_polygon_from_high_order_quad(e, polygon);
//                 // make_polygon_from_quad4(e, polygon);
//                 // } else {
//                 // 	make_polygon_from_curved_quad8(e, polygon);
//                 // }
//                 break;
//             }

//             default: {
//                 assert(false);
//                 break;
//             }
//         }
//     }

//     void make_polygon_3(const libMesh::Elem &e, libMesh::DenseMatrix<libMesh::Real> &polygon) {
//         auto n_nodes = e.n_nodes();

//         if (e.has_affine_map()) {
//             if (is_tri(e.type())) {
//                 n_nodes = 3;
//             } else if (is_quad(e.type())) {
//                 n_nodes = 4;
//             } else {
//                 assert(false && "handle special case");
//             }
//         }

//         polygon.resize(n_nodes, 3);

//         for (int i = 0; i < n_nodes; ++i) {
//             for (int j = 0; j < 3; ++j) {
//                 polygon(i, j) = e.point(i)(j);
//             }
//         }
//     }

//     // template<class FE>
//     // void mortar_assemble_weights_aux(const FE &fe, libMesh::DenseMatrix<libMesh::Real> &weights)
//     // {
//     // 	libMesh::DenseMatrix<libMesh::Real> elmat;
//     // 	elmat.resize(fe.get_phi().size(), fe.get_phi().size());
//     // 	elmat.zero();

//     // 	weights.resize(elmat.m(), elmat.n());
//     // 	weights.zero();

//     // 	const auto &test = fe.get_phi();
//     // 	const auto &JxW   = fe.get_JxW();

//     // 	const uint n_test  = test.size();
//     // 	const uint n_qp    = test[0].size();

//     // 	std::cout << n_qp << std::endl;

//     // 	for(uint qp = 0; qp < n_qp; ++qp) {
//     // 		for(uint i = 0; i < n_test; ++i) {
//     // 			for(uint j = 0; j < n_test; ++j) {
//     // 				elmat(i, j) += contract(test[i][qp], test[j][qp]) * JxW[qp];
//     // 			}
//     // 		}
//     // 	}

//     // 	libMesh::DenseVector<libMesh::Real> sum_elmat(n_test);
//     // 	sum_elmat.zero();
//     // 	libMesh::DenseVector<libMesh::Real> rhs(n_test);
//     // 	rhs.zero();

//     // 	libMesh::DenseVector<libMesh::Real> sol(n_test);
//     // 	sol.zero();

//     // 	for(uint i = 0; i < n_test; ++i) {
//     // 		for(uint j = 0; j < n_test; ++j) {
//     // 			sum_elmat(i) += elmat(i, j);
//     // 		}

//     // 		if(std::abs(sum_elmat(i)) < 1e-16) {
//     // 			sum_elmat(i) = 0;
//     // 			//set identity row where not defined
//     // 			for(uint j = 0; j < n_test; ++j) {
//     // 				elmat(i, j) = (i == j);
//     // 			}
//     // 		}
//     // 	}

//     // 	// std::cout << "-----------------------\n";
//     // 	// std::cout << "-----------------------\n";

//     // 	// elmat.print(std::cout);

//     // 	// std::cout << "-----------------------\n";

//     // 	for(uint i = 0; i < n_test; ++i) {
//     // 		if(sum_elmat(i) == 0) {
//     // 			continue;
//     // 		}

//     // 		rhs(i) = sum_elmat(i);

//     // 		elmat.cholesky_solve(rhs, sol);

//     // 		for(uint j = 0; j < n_test; ++j) {
//     // 			weights(i, j) = sol(j);
//     // 		}

//     // 		rhs(i) = 0;
//     // 	}

//     // 	//normalization for consistently scaled coefficients
//     // 	for(uint i = 0; i < n_test; ++i) {
//     // 		if(sum_elmat(i) == 0) {
//     // 			continue;
//     // 		}

//     // 		libMesh::Real t = 0;
//     // 		for(uint j = 0; j < n_test; ++j) {
//     // 			t += weights(i, j);
//     // 		}

//     // 		for(uint j = 0; j < n_test; ++j) {
//     // 			weights(i, j) *= 1./t;
//     // 		}
//     // 	}

//     // 	weights.print(std::cout);
//     // }

//     static libMesh::Real sum(const libMesh::Real &val) { return val; }

//     static libMesh::Real sum(const libMesh::VectorValue<libMesh::Real> &val) { return val(0) + val(1) + val(2); }

//     template <class FE>
//     void mortar_assemble_weights_aux(const FE &fe, libMesh::DenseMatrix<libMesh::Real> &weights) {
//         const auto &test = fe.get_phi();
//         const auto &JxW = fe.get_JxW();

//         const uint n_test = test.size();
//         const uint n_qp = test[0].size();

//         libMesh::DenseMatrix<libMesh::Real> elmat;
//         elmat.resize(n_test, n_test);
//         elmat.zero();

//         libMesh::DenseVector<libMesh::Real> sum_elmat(n_test);
//         sum_elmat.zero();

//         weights.resize(elmat.m(), elmat.n());
//         weights.zero();

//         for (uint qp = 0; qp < n_qp; ++qp) {
//             for (uint i = 0; i < n_test; ++i) {
//                 const auto val = sum(test[i][qp]) * JxW[qp];
//                 sum_elmat(i) += val;

//                 for (uint j = 0; j < n_test; ++j) {
//                     elmat(i, j) += contract(test[i][qp], test[j][qp]) * JxW[qp];
//                 }
//             }
//         }

//         libMesh::DenseVector<libMesh::Real> rhs(n_test);
//         rhs.zero();

//         libMesh::DenseVector<libMesh::Real> sol(n_test);
//         sol.zero();

//         for (uint i = 0; i < n_test; ++i) {
//             // for(uint j = 0; j < n_test; ++j) {
//             // 	sum_elmat(i) += elmat(i, j);
//             // }

//             if (std::abs(sum_elmat(i)) < 1e-16) {
//                 sum_elmat(i) = 0;
//                 // set identity row where not defined
//                 for (uint j = 0; j < n_test; ++j) {
//                     elmat(i, j) = (i == j);
//                 }
//             }
//         }

//         // std::cout << "-----------------------\n";

//         // std::cout << n_qp << std::endl;
//         // std::cout << "-----------------------\n";

//         // elmat.print(std::cout);

//         // std::cout << "-----------------------\n";

//         for (uint i = 0; i < n_test; ++i) {
//             if (sum_elmat(i) == 0) {
//                 continue;
//             }

//             rhs(i) = sum_elmat(i);

//             elmat.cholesky_solve(rhs, sol);

//             for (uint j = 0; j < n_test; ++j) {
//                 weights(i, j) = sol(j);
//             }

//             rhs(i) = 0;
//         }

//         // normalization for consistently scaled coefficients
//         for (uint i = 0; i < n_test; ++i) {
//             if (sum_elmat(i) == 0) {
//                 continue;
//             }

//             libMesh::Real t = 0;
//             for (uint j = 0; j < n_test; ++j) {
//                 t += weights(i, j);
//             }

//             for (uint j = 0; j < n_test; ++j) {
//                 weights(i, j) *= 1. / t;
//             }
//         }

//         // sum_elmat.print(std::cout);

//         // std::cout << "-----------------------\n";
//         // weights.print(std::cout);
//     }

//     void mortar_assemble_weights(const libMesh::FEVectorBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights) {
//         mortar_assemble_weights_aux(fe, weights);
//     }

//     void mortar_assemble_weights(const libMesh::FEBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights) {
//         mortar_assemble_weights_aux(fe, weights);
//     }

//     template <class FE>
//     void mortar_assemble_weighted_aux(const FE &trial_fe,
//                                       const FE &test_fe,
//                                       const libMesh::DenseMatrix<libMesh::Real> &weights,
//                                       libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
//             elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
//             elmat.zero();
//         }

//         const auto &trial = trial_fe.get_phi();
//         const auto &test = test_fe.get_phi();
//         const auto &JxW = test_fe.get_JxW();

//         const uint n_test = test.size();
//         const uint n_trial = trial.size();
//         const uint n_qp = test[0].size();

//         for (uint i = 0; i < n_test; ++i) {
//             for (uint qp = 0; qp < n_qp; ++qp) {
//                 auto w_test = test[i][qp] * 0.;

//                 for (uint k = 0; k < n_test; ++k) {
//                     w_test += test[k][qp] * weights(i, k);
//                 }

//                 for (uint j = 0; j < n_trial; ++j) {
//                     // assert(  JxW[qp] >= 0.0 );
//                     elmat(i, j) += contract(w_test, trial[j][qp]) * JxW[qp];
//                 }
//             }
//         }
//     }

//     void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
//                                          const libMesh::FEBase &test_fe,
//                                          const libMesh::DenseMatrix<libMesh::Real> &weights,
//                                          libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         mortar_assemble_weighted_aux(trial_fe, test_fe, weights, elmat);
//     }

//     void mortar_assemble_weighted_biorth(const libMesh::FEVectorBase &trial_fe,
//                                          const libMesh::FEVectorBase &test_fe,
//                                          const libMesh::DenseMatrix<libMesh::Real> &weights,
//                                          libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         mortar_assemble_weighted_aux(trial_fe, test_fe, weights, elmat);
//     }

//     // void integrate_scalar_function(
//     //     const libMesh::FEBase &test_fe,
//     //     const std::vector<double> &fun,
//     //     libMesh::DenseVector<libMesh::Real> &result
//     // )
//     // {
//     //     const auto &phi = test_fe.get_phi();
//     //     const auto &dx = test_fe.get_JxW();
//     //     const auto n_qp = fun.size();
//     //     const auto n_shape_functions = phi.size();

//     //     assert(n_qp == phi[0].size());
//     //     assert(n_qp == dx.size());

//     //     result.resize(n_shape_functions);
//     //     result.zero();

//     //     for(std::size_t i = 0; i < n_shape_functions; ++i) {
//     //         for(std::size_t qp = 0; qp < n_qp; ++qp) {
//     //             result(i) += phi[i][qp] * fun[qp] * dx[qp];
//     //         }
//     //     }
//     // }

//     void integrate_point_function(const int dim,
//                                   const libMesh::FEBase &test_fe,
//                                   const std::vector<libMesh::Point> &fun,
//                                   libMesh::DenseMatrix<libMesh::Real> &result) {
//         const auto &phi = test_fe.get_phi();
//         const auto &dx = test_fe.get_JxW();
//         const auto n_qp = fun.size();
//         const auto n_shape_functions = phi.size();

//         assert(n_qp == phi[0].size());
//         assert(n_qp == dx.size());

//         result.resize(n_shape_functions, dim);
//         result.zero();

//         for (std::size_t i = 0; i < n_shape_functions; ++i) {
//             for (std::size_t qp = 0; qp < n_qp; ++qp) {
//                 for (std::size_t d = 0; d < dim; ++d) {
//                     result(i, d) += phi[i][qp] * fun[qp](d) * dx[qp];
//                 }
//             }
//         }
//     }

//     void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
//                                          const libMesh::DenseMatrix<libMesh::Real> &trafo,
//                                          const libMesh::FEBase &test_fe,
//                                          const libMesh::DenseMatrix<libMesh::Real> &weights,
//                                          libMesh::DenseMatrix<libMesh::Real> &elmat) {
//         if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
//             elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
//             elmat.zero();
//         }

//         const auto &trial = trial_fe.get_phi();
//         const auto &test = test_fe.get_phi();
//         const auto &JxW = test_fe.get_JxW();

//         const uint n_test = test.size();
//         const uint n_trial = trial.size();
//         const uint n_qp = test[0].size();

//         std::vector<double> w_trial(n_trial);
//         std::vector<double> w_test(n_test);

//         for (uint qp = 0; qp < n_qp; ++qp) {
//             // build trial test
//             for (uint j = 0; j < n_trial; ++j) {
//                 auto w_trial_j = 0.;

//                 for (uint k = 0; k < n_trial; ++k) {
//                     w_trial_j += trial[k][qp] * trafo(j, k);
//                 }

//                 w_trial[j] = w_trial_j;
//             }

//             // build weigthed test
//             for (uint j = 0; j < n_test; ++j) {
//                 auto w_test_j = 0.;

//                 for (uint k = 0; k < n_test; ++k) {
//                     w_test_j += test[k][qp] * weights(j, k);
//                 }

//                 w_test[j] = w_test_j;
//             }

//             for (uint i = 0; i < n_test; ++i) {
//                 for (uint j = 0; j < n_trial; ++j) {
//                     elmat(i, j) += w_test[i] * w_trial[j] * JxW[qp];
//                 }
//             }
//         }
//     }
// }  // namespace utopia
