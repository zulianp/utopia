#include <libmesh/fe.h>
#include "MortarAssemble.hpp"
#include "utopia_Polygon.hpp"
#include "utopia_libmesh_Transform.hpp"
#include "utopia_libmesh_Utils.hpp"
#include "utopia_triangulate.hpp"

#include "Box.hpp"

#include <assert.h>
#include <algorithm>
#include <memory>
#include <numeric>

#include "utopia_intersector.hpp"

namespace utopia {
    template <class Left, class Right>
    libMesh::Real contract(const Left &left, const Right &right) {
        return left.contract(right);
    }

    libMesh::Real contract(const libMesh::Real &left, const libMesh::Real &right) { return left * right; }

    void compute_side_normal(const int dim, const libMesh::Elem &side, libMesh::Point &n) {
        using namespace libMesh;
        Point o, u, v;

        if (dim == 2) {
            assert(side.n_nodes() >= 2);
            o = side.point(0);
            u = side.point(1);
            u -= o;
            n(0) = u(1);
            n(1) = -u(0);

        } else {
            assert(dim >= 3);
            o = side.point(0);
            u = side.point(1);
            v = side.point(2);
            u -= o;
            v -= o;
            n = u.cross(v);
        }

        n *= 1. / n.norm();
    }

    int order_for_l2_integral(const int dim,
                              const libMesh::Elem &master_el,
                              const int master_order,
                              const libMesh::Elem &slave_el,
                              const int slave_order) {
        bool slave_has_affine_map = slave_el.has_affine_map();

        int order = 0;
        if (dim == 2 || dim == 1) {
            const auto actual_slave_order = slave_order * (is_quad(slave_el.type()) ? 2 : 1);
            order = master_order * (is_quad(master_el.type()) ? 2 : 1) + actual_slave_order;
            if (slave_has_affine_map) {
                order += std::max(actual_slave_order - 1, 0) * 2;
            }

        } else if (dim == 3) {
            const auto actual_slave_order = slave_order * (is_hex(slave_el.type()) ? 3 : 1);
            order = master_order * (is_hex(master_el.type()) ? 3 : 1) + actual_slave_order;
            if (slave_has_affine_map) {
                order += (actual_slave_order - 1) * 2;
            }
        } else {
            assert(false && "not supported yet for dim != 2 or dim != 3");
        }

        return order;
    }

    void print(const libMesh::QBase &ir, std::ostream &os) {
        os << "points:\n[\n";
        for (int i = 0; i < ir.n_points(); ++i) {
            os << "\t" << ir.qp(i)(0) << ", " << ir.qp(i)(1) << ", " << ir.qp(i)(2) << "\n";
        }

        os << "]\n";

        os << "weights:\n[\n\t";
        for (int i = 0; i < ir.n_points(); ++i) {
            os << ir.w(i) << " ";
        }
        os << "\n]\n";
    }

    double sum_of_weights(const libMesh::QBase &ir) {
        double ret = 0;
        for (int i = 0; i < ir.n_points(); ++i) {
            ret += ir.w(i);
        }
        return ret;
    }

    template <class V, class T>
    static void add(const V &p, const T &alpha, const V &v, V &result) {
        assert(result.size() == p.size());

        for (int i = 0; i < result.size(); ++i) {
            result(i) = p(i) + alpha * v(i);
        }
    }

    // bool biorthgonal_weights(const int type, libMesh::Real &w_ii, libMesh::Real &w_ij) {
    //     using namespace libMesh;

    //     switch (type) {
    //         case EDGE2: {
    //             w_ii = 2.0;
    //             w_ij = -1.0;
    //             return true;
    //         }

    //         case TRI3: {
    //             w_ii = 3.0;
    //             w_ij = -1.0;
    //             return true;
    //         }

    //         case TET4: {
    //             w_ii = 4.0;
    //             w_ij = -1.0;
    //             return true;
    //         }
    //             // These do not work:
    //             // case QUAD4:
    //             // {
    //             //  w_ii = 4.0;
    //             //  w_ij = -1.0;
    //             //  return true;
    //             // }

    //             // case HEX8:
    //             // {
    //             //  w_ii = 8.0;
    //             //  w_ij = -1.0;
    //             //  return true;
    //             // }

    //         default: {
    //             w_ii = 1.0;
    //             w_ij = 0.0;

    //             static bool error_msg_printed = false;

    //             if (!error_msg_printed) {
    //                 std::cerr << "[Error] biorthgonal weights not supported for element type: " << type << std::endl;
    //                 error_msg_printed = true;
    //             }

    //             assert(false && "TODO: add the weights for the missing element");
    //             return false;
    //         }
    //     }
    // }

    // template<class FE>
    // void mortar_assemble_weights_aux(const FE &fe, libMesh::DenseMatrix<libMesh::Real> &weights)
    // {
    //  libMesh::DenseMatrix<libMesh::Real> elmat;
    //  elmat.resize(fe.get_phi().size(), fe.get_phi().size());
    //  elmat.zero();

    //  weights.resize(elmat.m(), elmat.n());
    //  weights.zero();

    //  const auto &test = fe.get_phi();
    //  const auto &JxW   = fe.get_JxW();

    //  const uint n_test  = test.size();
    //  const uint n_qp    = test[0].size();

    //  std::cout << n_qp << std::endl;

    //  for(uint qp = 0; qp < n_qp; ++qp) {
    //      for(uint i = 0; i < n_test; ++i) {
    //          for(uint j = 0; j < n_test; ++j) {
    //              elmat(i, j) += contract(test[i][qp], test[j][qp]) * JxW[qp];
    //          }
    //      }
    //  }

    //  libMesh::DenseVector<libMesh::Real> sum_elmat(n_test);
    //  sum_elmat.zero();
    //  libMesh::DenseVector<libMesh::Real> rhs(n_test);
    //  rhs.zero();

    //  libMesh::DenseVector<libMesh::Real> sol(n_test);
    //  sol.zero();

    //  for(uint i = 0; i < n_test; ++i) {
    //      for(uint j = 0; j < n_test; ++j) {
    //          sum_elmat(i) += elmat(i, j);
    //      }

    //      if(std::abs(sum_elmat(i)) < 1e-16) {
    //          sum_elmat(i) = 0;
    //          //set identity row where not defined
    //          for(uint j = 0; j < n_test; ++j) {
    //              elmat(i, j) = (i == j);
    //          }
    //      }
    //  }

    //  // std::cout << "-----------------------\n";
    //  // std::cout << "-----------------------\n";

    //  // elmat.print(std::cout);

    //  // std::cout << "-----------------------\n";

    //  for(uint i = 0; i < n_test; ++i) {
    //      if(sum_elmat(i) == 0) {
    //          continue;
    //      }

    //      rhs(i) = sum_elmat(i);

    //      elmat.cholesky_solve(rhs, sol);

    //      for(uint j = 0; j < n_test; ++j) {
    //          weights(i, j) = sol(j);
    //      }

    //      rhs(i) = 0;
    //  }

    //  //normalization for consistently scaled coefficients
    //  for(uint i = 0; i < n_test; ++i) {
    //      if(sum_elmat(i) == 0) {
    //          continue;
    //      }

    //      libMesh::Real t = 0;
    //      for(uint j = 0; j < n_test; ++j) {
    //          t += weights(i, j);
    //      }

    //      for(uint j = 0; j < n_test; ++j) {
    //          weights(i, j) *= 1./t;
    //      }
    //  }

    //  weights.print(std::cout);
    // }

    static libMesh::Real sum(const libMesh::Real &val) { return val; }

    static libMesh::Real sum(const libMesh::VectorValue<libMesh::Real> &val) { return val(0) + val(1) + val(2); }

    template <class FE>
    void mortar_assemble_weights_aux(const FE &fe, libMesh::DenseMatrix<libMesh::Real> &weights) {
        const auto &test = fe.get_phi();
        const auto &JxW = fe.get_JxW();

        const uint n_test = test.size();
        const uint n_qp = test[0].size();

        libMesh::DenseMatrix<libMesh::Real> elmat;
        elmat.resize(n_test, n_test);
        elmat.zero();

        libMesh::DenseVector<libMesh::Real> sum_elmat(n_test);
        sum_elmat.zero();

        weights.resize(elmat.m(), elmat.n());
        weights.zero();

        for (uint qp = 0; qp < n_qp; ++qp) {
            for (uint i = 0; i < n_test; ++i) {
                const auto val = sum(test[i][qp]) * JxW[qp];
                sum_elmat(i) += val;

                for (uint j = 0; j < n_test; ++j) {
                    elmat(i, j) += contract(test[i][qp], test[j][qp]) * JxW[qp];
                }
            }
        }

        libMesh::DenseVector<libMesh::Real> rhs(n_test);
        rhs.zero();

        libMesh::DenseVector<libMesh::Real> sol(n_test);
        sol.zero();

        for (uint i = 0; i < n_test; ++i) {
            // for(uint j = 0; j < n_test; ++j) {
            //  sum_elmat(i) += elmat(i, j);
            // }

            if (std::abs(sum_elmat(i)) < 1e-16) {
                sum_elmat(i) = 0;
                // set identity row where not defined
                for (uint j = 0; j < n_test; ++j) {
                    elmat(i, j) = (i == j);
                }
            }
        }

        // std::cout << "-----------------------\n";

        // std::cout << n_qp << std::endl;
        // std::cout << "-----------------------\n";

        // elmat.print(std::cout);

        // std::cout << "-----------------------\n";

        for (uint i = 0; i < n_test; ++i) {
            if (sum_elmat(i) == 0) {
                continue;
            }

            rhs(i) = sum_elmat(i);

            elmat.cholesky_solve(rhs, sol);

            for (uint j = 0; j < n_test; ++j) {
                weights(i, j) = sol(j);
            }

            rhs(i) = 0;
        }

        // normalization for consistently scaled coefficients
        for (uint i = 0; i < n_test; ++i) {
            if (sum_elmat(i) == 0) {
                continue;
            }

            libMesh::Real t = 0;
            for (uint j = 0; j < n_test; ++j) {
                t += weights(i, j);
            }

            for (uint j = 0; j < n_test; ++j) {
                weights(i, j) *= 1. / t;
            }
        }

        // sum_elmat.print(std::cout);

        // std::cout << "-----------------------\n";
        // weights.print(std::cout);
    }

    void mortar_assemble_weights(const libMesh::FEVectorBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights) {
        mortar_assemble_weights_aux(fe, weights);
    }

    void mortar_assemble_weights(const libMesh::FEBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights) {
        mortar_assemble_weights_aux(fe, weights);
    }

    template <class FE>
    void mortar_assemble_weighted_aux(const FE &trial_fe,
                                      const FE &test_fe,
                                      const libMesh::DenseMatrix<libMesh::Real> &weights,
                                      libMesh::DenseMatrix<libMesh::Real> &elmat) {
        if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
            elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
            elmat.zero();
        }

        const auto &trial = trial_fe.get_phi();
        const auto &test = test_fe.get_phi();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_trial = trial.size();
        const uint n_qp = test[0].size();

        for (uint i = 0; i < n_test; ++i) {
            for (uint qp = 0; qp < n_qp; ++qp) {
                auto w_test = test[i][qp] * 0.;

                for (uint k = 0; k < n_test; ++k) {
                    w_test += test[k][qp] * weights(i, k);
                }

                for (uint j = 0; j < n_trial; ++j) {
                    // assert(  JxW[qp] >= 0.0 );
                    elmat(i, j) += contract(w_test, trial[j][qp]) * JxW[qp];
                }
            }
        }
    }

    void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
                                         const libMesh::FEBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat) {
        mortar_assemble_weighted_aux(trial_fe, test_fe, weights, elmat);
    }

    void mortar_assemble_weighted_biorth(const libMesh::FEVectorBase &trial_fe,
                                         const libMesh::FEVectorBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat) {
        mortar_assemble_weighted_aux(trial_fe, test_fe, weights, elmat);
    }

    // void integrate_scalar_function(
    //     const libMesh::FEBase &test_fe,
    //     const std::vector<double> &fun,
    //     libMesh::DenseVector<libMesh::Real> &result
    // )
    // {
    //     const auto &phi = test_fe.get_phi();
    //     const auto &dx = test_fe.get_JxW();
    //     const auto n_qp = fun.size();
    //     const auto n_shape_functions = phi.size();

    //     assert(n_qp == phi[0].size());
    //     assert(n_qp == dx.size());

    //     result.resize(n_shape_functions);
    //     result.zero();

    //     for(std::size_t i = 0; i < n_shape_functions; ++i) {
    //         for(std::size_t qp = 0; qp < n_qp; ++qp) {
    //             result(i) += phi[i][qp] * fun[qp] * dx[qp];
    //         }
    //     }
    // }

    void integrate_point_function(const int dim,
                                  const libMesh::FEBase &test_fe,
                                  const std::vector<libMesh::Point> &fun,
                                  libMesh::DenseMatrix<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = fun.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());
        assert(n_qp == dx.size());

        result.resize(n_shape_functions, dim);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                for (std::size_t d = 0; d < dim; ++d) {
                    result(i, d) += phi[i][qp] * fun[qp](d) * dx[qp];
                }
            }
        }
    }

    // void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEVectorBase &test_fe,
    //                                                     const int dim,
    //                                                     const libMesh::Point &surf_normal,
    //                                                     const libMesh::Point &plane_normal,
    //                                                     const libMesh::Real &plane_offset,
    //                                                     const libMesh::DenseMatrix<libMesh::Real> &weights,
    //                                                     libMesh::DenseMatrix<libMesh::Real> &normals,
    //                                                     libMesh::DenseVector<libMesh::Real> &gap) {
    //     using namespace libMesh;

    //     if (normals.m() != test_fe.get_phi().size() / dim || dim != normals.n()) {
    //         normals.resize(test_fe.get_phi().size() / dim, dim);
    //         gap.resize(test_fe.get_phi().size());
    //     }

    //     normals.zero();
    //     gap.zero();

    //     const auto &test = test_fe.get_phi();
    //     const auto &point = test_fe.get_xyz();
    //     const auto &JxW = test_fe.get_JxW();

    //     const uint n_test = test.size();
    //     const uint n_qp = test[0].size();

    //     DenseVector<Real> p(dim);
    //     DenseVector<Real> v(dim);

    //     DenseVector<Real> s_n(dim);
    //     DenseVector<Real> p_n(dim);

    //     for (uint i = 0; i < dim; ++i) {
    //         p_n(i) = plane_normal(i);
    //         s_n(i) = surf_normal(i);
    //     }

    //     for (uint qp = 0; qp < n_qp; ++qp) {
    //         p(0) = point[qp](0);
    //         p(1) = point[qp](1);

    //         if (dim > 2) {
    //             p(2) = point[qp](2);
    //         }

    //         Real isect = 0;

    //         Intersector::intersect_ray_with_plane(
    //             dim, 1, &p.get_values()[0], &s_n.get_values()[0], &p_n.get_values()[0], plane_offset, &isect);

    //         v = s_n;
    //         v *= isect;
    //         // quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]);

    //         for (uint i = 0; i < n_test; ++i) {
    //             auto biorth_test = weights(i, 0) * test[0][qp];

    //             for (uint k = 0; k < test.size(); ++k) {
    //                 biorth_test += weights(i, k) * test[k][qp];
    //             }

    //             gap(i) += biorth_test(0) * isect * JxW[qp];

    //             for (uint d = 0; d < dim; ++d) {
    //                 normals.get_values()[i] += biorth_test(d) * surf_normal(d) * JxW[qp];
    //             }
    //         }
    //     }

    //     // gap.print(std::cout);
    // }

    template <typename T>
    void convert_point_to_vector(const int dim, const libMesh::Point &point, std::vector<T> &point_vec) {
        point_vec.resize(dim);
        for (int i = 0; i < dim; ++i) {
            point_vec[i] = point(i);
        }
    }

    // void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEBase &test_fe,
    //                                                     const int dim,
    //                                                     const libMesh::Point &surf_normal,
    //                                                     const libMesh::Point &plane_normal,
    //                                                     const libMesh::Real &plane_offset,
    //                                                     const libMesh::DenseMatrix<libMesh::Real> &weights,
    //                                                     libMesh::DenseMatrix<libMesh::Real> &normals,
    //                                                     libMesh::DenseVector<libMesh::Real> &gap,
    //                                                     const bool visdbg) {
    //     using namespace libMesh;

    //     if (normals.m() != test_fe.get_phi().size() || dim != normals.n()) {
    //         normals.resize(test_fe.get_phi().size(), dim);
    //         normals.zero();
    //         gap.resize(test_fe.get_phi().size());
    //         gap.zero();
    //     }

    //     const auto &test = test_fe.get_phi();
    //     // const auto &grad   = test_fe.get_dphi();
    //     const auto &point = test_fe.get_xyz();
    //     const auto &JxW = test_fe.get_JxW();

    //     const uint n_test = test.size();
    //     const uint n_qp = test[0].size();

    //     std::vector<Real> surf_normal_v, plane_normal_v;
    //     convert_point_to_vector(dim, surf_normal, surf_normal_v);
    //     convert_point_to_vector(dim, plane_normal, plane_normal_v);

    //     DenseVector<Real> p(dim);
    //     // DenseMatrix<Real> v(dim); //visdbg

    //     for (uint i = 0; i < n_test; ++i) {
    //         for (uint qp = 0; qp < n_qp; ++qp) {
    //             p(0) = point[qp](0);
    //             p(1) = point[qp](1);

    //             if (dim > 2) {
    //                 p(2) = point[qp](2);
    //             }

    //             Real isect = 0;
    //             Intersector::intersect_ray_with_plane(
    //                 dim, 1, &p.get_values()[0], &surf_normal_v[0], &plane_normal_v[0], plane_offset, &isect);
    //             // assert(isect > 0);
    //             // if(visdbg) {
    //             //  v.get_values() = surf_normal_v; //visdbg
    //             //  v *= isect; //visdbg
    //             //  quiver(dim, 1, &p.get_values()[0], &v.get_values()[0]); //visdbg
    //             // }

    //             auto biorth_test = weights(i, 0) * test[0][qp];

    //             for (uint k = 1; k < n_test; ++k) {
    //                 biorth_test += weights(i, k) * test[k][qp];
    //             }

    //             gap(i) += biorth_test * isect * JxW[qp];

    //             for (uint d = 0; d < dim; ++d) {
    //                 normals(i, d) += biorth_test * surf_normal(d) * JxW[qp];
    //             }
    //         }
    //     }
    // }

    void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &trafo,
                                         const libMesh::FEBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat) {
        if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
            elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
            elmat.zero();
        }

        const auto &trial = trial_fe.get_phi();
        const auto &test = test_fe.get_phi();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_trial = trial.size();
        const uint n_qp = test[0].size();

        std::vector<double> w_trial(n_trial);
        std::vector<double> w_test(n_test);

        for (uint qp = 0; qp < n_qp; ++qp) {
            // build trial test
            for (uint j = 0; j < n_trial; ++j) {
                auto w_trial_j = 0.;

                for (uint k = 0; k < n_trial; ++k) {
                    w_trial_j += trial[k][qp] * trafo(j, k);
                }

                w_trial[j] = w_trial_j;
            }

            // build weigthed test
            for (uint j = 0; j < n_test; ++j) {
                auto w_test_j = 0.;

                for (uint k = 0; k < n_test; ++k) {
                    w_test_j += test[k][qp] * weights(j, k);
                }

                w_test[j] = w_test_j;
            }

            for (uint i = 0; i < n_test; ++i) {
                for (uint j = 0; j < n_trial; ++j) {
                    elmat(i, j) += w_test[i] * w_trial[j] * JxW[qp];
                }
            }
        }
    }

    double ref_volume(int type) {
        if (is_hex(type)) {
            return 8.;
        } else if (is_tet(type)) {
            return 1. / 6.;
        } else if (is_quad(type)) {
            return 4.;
        } else if (is_edge(type)) {
            return 2.;
        } else if (is_tri(type)) {
            return 0.5;
        } else if (is_prism(type)) {
            return 1.;
        } else if (is_pyramid(type)) {
            return 1. / 0.75;
        } else {
            assert(false);
            return 1.;
        }
    }

    double ref_area_of_surf(int type) {
        if (is_tri(type)) {
            return 2.;
        } else if (is_quad(type)) {
            return 2.;
        } else if (is_hex(type)) {
            return 1.;
        } else if (is_tet(type)) {
            return 0.5;
        } else if (is_edge(type)) {
            // Is this correct?
            assert(false);
            return 0.5;
        } else {
            assert(false && "add special case");
            return 1.;
        }
    }

    template <class FE>
    void mortar_assemble_aux(const FE &trial_fe, const FE &test_fe, libMesh::DenseMatrix<libMesh::Real> &elmat) {
        if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
            elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
            elmat.zero();
        }

        const auto &trial = trial_fe.get_phi();
        const auto &test = test_fe.get_phi();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_trial = trial.size();
        const uint n_qp = test[0].size();

        assert(test[0].size() == trial[0].size());

        for (uint i = 0; i < n_test; ++i) {
            for (uint j = 0; j < n_trial; ++j) {
                for (uint qp = 0; qp < n_qp; ++qp) {
                    elmat(i, j) += contract(test[i][qp], trial[j][qp]) * JxW[qp];
                }
            }
        }
    }

    libMesh::Real len(const libMesh::Real val) { return std::abs(val); }

    template <class Vec>
    libMesh::Real len(const Vec &val) {
        return val.size();
    }

    static inline bool is_vec(const libMesh::Real val) { return false; }

    template <class Vec>
    bool is_vec(const Vec &val) {
        return true;
    }

    // void make_tp(const int i, libMesh::Real &val) {}

    // template <class Vec>
    // void make_tp(const int i, Vec &val) {
    //     libMesh::Real s = val(i);
    //     val.zero();
    //     val(i) = s;
    // }

    template <class FE>
    void mortar_assemble_biorth_aux(const FE &trial_fe,
                                    const FE &test_fe,
                                    const libMesh::Real &wii,
                                    const libMesh::Real &wij,
                                    libMesh::DenseMatrix<libMesh::Real> &elmat) {
        if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
            elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
            elmat.zero();
        }

        const auto &trial = trial_fe.get_phi();
        const auto &test = test_fe.get_phi();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_trial = trial.size();
        const uint n_qp = test[0].size();

        for (uint qp = 0; qp < n_qp; ++qp) {
            for (uint i = 0; i < n_test; ++i) {
                auto biorth_test = ((0 == i) ? wii : wij) * test[0][qp];

                for (uint k = 1; k < n_test; ++k) {
                    biorth_test += ((k == i) ? wii : wij) * test[k][qp];
                }

                for (uint j = 0; j < n_trial; ++j) {
                    elmat(i, j) += contract(biorth_test, trial[j][qp]) * JxW[qp];
                }
            }
        }
    }

    void mortar_assemble(const libMesh::FEBase &trial_fe,
                         const libMesh::FEBase &test_fe,
                         libMesh::DenseMatrix<libMesh::Real> &elmat) {
        mortar_assemble_aux(trial_fe, test_fe, elmat);
    }

    void mortar_assemble(const libMesh::FEVectorBase &trial_fe,
                         const libMesh::FEVectorBase &test_fe,
                         libMesh::DenseMatrix<libMesh::Real> &elmat) {
        mortar_assemble_aux(trial_fe, test_fe, elmat);
    }

    void mortar_assemble_biorth(const libMesh::FEBase &trial_fe,
                                const libMesh::FEBase &test_fe,
                                const int type,
                                libMesh::DenseMatrix<libMesh::Real> &elmat) {
        libMesh::Real w_ii, w_ij;
        biorthgonal_weights(type, w_ii, w_ij);
        mortar_assemble_biorth_aux(trial_fe, test_fe, w_ii, w_ij, elmat);
    }

    void mortar_assemble_biorth(const libMesh::FEVectorBase &trial_fe,
                                const libMesh::FEVectorBase &test_fe,
                                const int type,
                                libMesh::DenseMatrix<libMesh::Real> &elmat) {
        libMesh::Real w_ii, w_ij;
        biorthgonal_weights(type, w_ii, w_ij);
        mortar_assemble_biorth_aux(trial_fe, test_fe, w_ii, w_ij, elmat);
    }

    template <class FE>
    void mortar_assemble_biorth_aux(const int dim,
                                    const FE &trial_fe,
                                    const FE &test_fe,
                                    const libMesh::Real &wii,
                                    const libMesh::Real &wij,
                                    const libMesh::DenseVector<libMesh::Real> &indicator,
                                    libMesh::DenseMatrix<libMesh::Real> &elmat) {
        if (elmat.m() != test_fe.get_phi().size() || elmat.n() != trial_fe.get_phi().size()) {
            elmat.resize(test_fe.get_phi().size(), trial_fe.get_phi().size());
            elmat.zero();
        }

        const auto &trial = trial_fe.get_phi();
        const auto &test = test_fe.get_phi();
        const auto &JxW = test_fe.get_JxW();

        const uint n_test = test.size();
        const uint n_trial = trial.size();
        const uint n_qp = test[0].size();

        //      bool v = is_vec(test[0][0]);

        for (uint qp = 0; qp < n_qp; ++qp) {
            for (uint i = 0; i < n_test; ++i) {
                auto biorth_test = ((0 == i) ? wii : wij) * test[0][qp];

                for (uint k = 1; k < n_test; ++k) {
                    biorth_test += ((k == i) ? wii : wij) * test[k][qp];
                }

                for (int k = 0; k < dim; ++k) {
                    if (i % dim == k) {
                        make_tp(k, biorth_test);
                        break;
                    }
                }

                // if(indicator(i) > 0)
                //  std::cout <<  biorth_test << std::endl;

                for (uint j = 0; j < n_trial; ++j) {
                    elmat(i, j) += indicator(i) * contract(biorth_test, trial[j][qp]) * JxW[qp];
                }
            }
        }
    }

    void mortar_assemble_biorth(const int dim,
                                const libMesh::FEBase &trial_fe,
                                const libMesh::FEBase &test_fe,
                                const int type,
                                const libMesh::DenseVector<libMesh::Real> &indicator,
                                libMesh::DenseMatrix<libMesh::Real> &elmat) {
        libMesh::Real w_ii, w_ij;
        biorthgonal_weights(type, w_ii, w_ij);
        mortar_assemble_biorth_aux(dim, trial_fe, test_fe, w_ii, w_ij, indicator, elmat);
    }

    void mortar_assemble_biorth(const int dim,
                                const libMesh::FEVectorBase &trial_fe,
                                const libMesh::FEVectorBase &test_fe,
                                const int type,
                                const libMesh::DenseVector<libMesh::Real> &indicator,
                                libMesh::DenseMatrix<libMesh::Real> &elmat) {
        libMesh::Real w_ii, w_ij;
        biorthgonal_weights(type, w_ii, w_ij);
        mortar_assemble_biorth_aux(dim, trial_fe, test_fe, w_ii, w_ij, indicator, elmat);
    }

    void transform_to_reference(const Transform &trans, const int type, const QMortar &global_ir, QMortar &ref_ir) {
        const int dim = global_ir.get_dim();
        libMesh::Point p;
        ref_ir.resize(global_ir.n_points());

        // libMesh::DenseMatrix<libMesh::Real> A_inv;

        const double factor = ref_volume(type);

        double sum_of_weights = 0.;
        for (int i = 0; i < global_ir.n_points(); ++i) {
            p(0) = global_ir.qp(i)(0);
            p(1) = global_ir.qp(i)(1);
            p(2) = global_ir.qp(i)(2);

            trans.transform_to_reference(p, ref_ir.get_points()[i]);

            ref_ir.get_weights()[i] = global_ir.w(i) * factor;

            // for debugging
            sum_of_weights += ref_ir.get_weights()[i];
        }

        // assert(sum_of_weights > 0.);
        // assert(sum_of_weights <= factor + 5e-4);
    }

    void transform_to_reference_surf(const Transform &trans,
                                     const int type,
                                     const QMortar &global_ir,
                                     QMortar &ref_ir) {
        assert(is_valid_elem_type(type));

        const int dim = global_ir.get_dim();
        libMesh::Point p;
        ref_ir.resize(global_ir.n_points());

        libMesh::DenseMatrix<libMesh::Real> A_inv;

        const double factor = ref_area_of_surf(type);

        for (int i = 0; i < global_ir.n_points(); ++i) {
            p(0) = global_ir.qp(i)(0);
            p(1) = global_ir.qp(i)(1);
            p(2) = global_ir.qp(i)(2);

            trans.transform_to_reference(p, ref_ir.get_points()[i]);
            ref_ir.get_weights()[i] = global_ir.w(i) * factor;
        }
    }

}  // namespace utopia
