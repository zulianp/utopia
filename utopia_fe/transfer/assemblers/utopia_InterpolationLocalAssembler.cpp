#include "utopia_InterpolationLocalAssembler.hpp"

namespace utopia {

    using Vector2 = Intersector::Vector2;
    inline static bool inside_half_plane(const Vector2 &e1, const Vector2 &e2, const Vector2 &point, const double tol)
    {
        const Vector2 u = e1 - e2;
        const Vector2 v = point - e2;

        const double dist = (u.x * v.y) - (v.x * u.y);
        return dist <= tol;
    }

    std::shared_ptr<Transform> InterpolationLocalAssembler::get_trafo(const Elem &elem) const
    {
        std::shared_ptr<Transform> elem_trafo;
        if(force_affine_ || elem.has_affine_map()) {
            if(dim == 1) {
                elem_trafo = std::make_shared<Transform1>(elem);
            } else  if(dim == 2) {
                if(elem.dim() == 1) {
                    elem_trafo = std::make_shared<Transform1>(elem);
                } else {
                    elem_trafo = std::make_shared<AffineTransform2>(elem);
                }
            } else {
                assert(dim == 3);

                if(elem.dim() == 2) {
                    elem_trafo = std::make_shared<Transform2>(elem);
                } else {
                    elem_trafo = std::make_shared<AffineTransform3>(elem);
                }
            }
        } else {
            if(dim == 1) {
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
        }

        return elem_trafo;
    }

    bool InterpolationLocalAssembler::assemble(
        const Elem &trial,
        FEType trial_type,
        const Elem &test,
        FEType test_type,
        Matrix &mat
        )
    {
        int n_potential_nodes = test.n_nodes();

        std::vector<int> test_dofs;
        contained_points(trial, test, test_dofs);

        if(test_dofs.empty()) return false;

        std::shared_ptr<Transform> trial_trafo = get_trafo(trial);
        std::shared_ptr<Transform> test_trafo  = get_trafo(test);

        init_q(test_dofs.size());

        int i_ref = 0;
        for(auto i : test_dofs) {
            trial_trafo->transform_to_reference(test.node_ref(i), q_trial->get_points()[i_ref++]);
        }

        i_ref = 0;
        for(auto i : test_dofs) {
            test_trafo->transform_to_reference(test.node_ref(i), q_test->get_points()[i_ref++]);
        }

        auto trial_fe = libMesh::FEBase::build(trial.dim(), trial_type);
        trial_fe->attach_quadrature_rule(q_trial.get());
        const auto &trial_shape_fun = trial_fe->get_phi();
        trial_fe->reinit(&trial);

        auto test_fe = libMesh::FEBase::build(test.dim(), test_type);
        test_fe->attach_quadrature_rule(q_test.get());
        const auto &test_shape_fun = test_fe->get_phi();
        test_fe->reinit(&test);

        mat.resize(test_shape_fun.size(), trial_shape_fun.size());

        for(std::size_t i = 0; i < test_shape_fun.size(); ++i) {
            for(std::size_t j = 0; j < trial_shape_fun.size(); ++j) {
                for(std::size_t k = 0; k < test_dofs.size(); ++k) {
                    auto tf = test_shape_fun.at(i).at(k);

                    if(tf < 0.9) {
                        //exploiting the lagrange property
                        tf = 0.;
                    }

                    mat(i, j) += trial_shape_fun.at(j).at(k) * tf;
                }
            }
        }

        // mat.print();

        assert((test_type != libMesh::FIRST || trial_type != libMesh::FIRST || check_valid(mat)));
        return true;
    }

    bool InterpolationLocalAssembler::check_valid(const Matrix &mat) const
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

    void InterpolationLocalAssembler::init_q(const std::size_t n_qp)
    {
        if(!q_trial) {
            q_trial = std::make_shared<QMortar>(dim);
            q_test  = std::make_shared<QMortar>(dim);
        }

        q_trial->resize(n_qp);
        q_test->resize(n_qp);

        std::fill(q_trial->get_weights().begin(), q_trial->get_weights().end(), 0.);
        std::fill(q_test->get_weights().begin(),  q_test->get_weights().end(), 0.);
    }

    void InterpolationLocalAssembler::contained_points(const Elem &trial, const Elem &test, std::vector<int> &test_dofs)
    {
        int n_potential_nodes = test.n_nodes();

        if(nested_meshes) {
            //check if there is an intersection then...


            test_dofs.resize(n_potential_nodes);

            for(int i = 0; i < n_potential_nodes; ++i) {
                test_dofs[i] = i;
            }

            return;
        }

        test_dofs.reserve(n_potential_nodes);

        if(dim == 1) {
            contained_points_1(trial, test, test_dofs);
            return;
        } else if(dim == 2) {
            contained_points_2(trial, test, test_dofs);
            return;
        } else if(dim == 3) {
            assert(dim == 3);
            contained_points_3(trial, test, test_dofs);
            return;
        }

        for(int i = 0; i < n_potential_nodes; ++i) {
            auto const & test_node = test.node_ref(i);
            if(trial.contains_point(test_node, tol)) {
                test_dofs.push_back(i);
            }
        }
    }

    void InterpolationLocalAssembler::contained_points_1(const Elem &trial, const Elem &test, std::vector<int> &test_dofs)
    {
        libMesh::Point p1 = trial.node_ref(0), p2 = trial.node_ref(1);
        libMesh::Point u = p2 - p1;
        auto len_u = u.norm();
        u /= len_u;

        test_dofs.clear();

        auto test_n_nodes = test.n_nodes();
        for(int i = 0; i < test_n_nodes; ++i) {
            const auto &q = test.node_ref(i);
            double d = u.contract(q - p1);
            if(d >= -tol && d <= len_u + tol) {
                test_dofs.push_back(i);
            }
        }
    }

    void InterpolationLocalAssembler::contained_points_2(const Elem &trial, const Elem &test, std::vector<int> &test_dofs)
    {
        make_polygon(trial, trial_pts);
        make_polygon(test,  test_pts);

        const auto &trial_poly = trial_pts.get_values();
        const int trial_n_nodes = trial_poly.size() / 2;
        const int test_n_nodes  = test.n_nodes();

        test_dofs.clear();

        Vector2 e1, e2, p, s, e;

        std::vector<bool> is_inside(test_n_nodes, true);

        for(int i = 0; i < trial_n_nodes; ++i) {
            const int i2x = i * 2;
            const int i2y = i2x + 1;

            const int i2p1x = 2 * (((i + 1) == trial_n_nodes)? 0 : (i + 1));
            const int i2p1y = i2p1x + 1;

            e1 = Vector2(trial_poly[i2x],   trial_poly[i2y]);
            e2 = Vector2(trial_poly[i2p1x], trial_poly[i2p1y]);

            for(int j = 0; j < test_n_nodes; ++j) {
                const auto &test_node = test.node_ref(j);
                p = Vector2(test_node(0), test_node(1));
                is_inside[j] = is_inside[j] && inside_half_plane(e1, e2, p, tol);
            }
        }

        for(int i = 0; i < test_n_nodes; ++i) {
            if(is_inside[i]) {
                test_dofs.push_back(i);
            }
        }
    }

    void InterpolationLocalAssembler::contained_points_3(const Elem &trial, const Elem &test, std::vector<int> &test_dofs) const
    {
        Polyhedron poly;
        make_polyhedron(trial, poly);

        auto n_half_spaces = poly.n_elements;

        std::vector<double> plane_normals(n_half_spaces * 3, 0.);
        std::vector<double> plane_dists_from_origin(n_half_spaces, 0.);
        Intersector::make_h_polyhedron_from_polyhedron(poly, &plane_normals[0], &plane_dists_from_origin[0]);

        double p[3];

        for(int k = 0; k < test.n_nodes(); ++k) {
            const auto & node = test.node_ref(k);
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
                test_dofs.push_back(k);
            }
        }
    }
}
