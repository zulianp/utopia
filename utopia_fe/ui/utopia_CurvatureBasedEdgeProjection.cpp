#include "utopia_CurvatureBasedEdgeProjection.hpp"
#include "utopia_NormalTangentialCoordinateSystem.hpp"
#include "utopia_LibMeshToMoonolithConvertions.hpp"

#include "moonolith_pn_triangle.hpp"

namespace utopia {

 static bool assemble_smooth_normals(const libMesh::MeshBase &mesh,
   const libMesh::DofMap &dof_map,
   const std::vector<int> &boundary_tags,
   UVector &normals)
 {
    using namespace libMesh;

    const uint spatial_dim = mesh.spatial_dimension();
    const uint n_dims = mesh.mesh_dimension();
    std::unique_ptr<FEBase> fe = FEBase::build(n_dims, dof_map.variable_order(0));
    fe->get_normals();
    fe->get_phi();
    fe->get_JxW();

    QGauss quad(n_dims-1 , FIFTH);


    const SizeType local_dofs = dof_map.n_local_dofs();
    normals = local_zeros(local_dofs * spatial_dim);

    std::vector<libMesh::dof_id_type> dof_indices;
    DenseVector<Real> vec, local_touched;

    SizeType n_detected_side_sets = 0;
        { //synch-block begin
            Write<UVector> w_n(normals);

            for(auto e_it = mesh.active_local_elements_begin(); e_it != mesh.active_local_elements_end(); ++e_it) {
                const auto &e = **e_it;

                for(uint side = 0; side < e.n_sides(); ++side) {
                    if(e.neighbor_ptr(side) != nullptr) {continue;}

                    //if empty take them all
                    bool select = boundary_tags.empty();
                    for(auto t : boundary_tags) {
                        if(mesh.get_boundary_info().has_boundary_id(&e, side, t)) {
                            select = true;
                            break;
                        }
                    }

                    if(!select) continue;
                    ++n_detected_side_sets;

                    fe->attach_quadrature_rule(&quad);
                    fe->reinit(&e, side);

                    const auto &fe_normals = fe->get_normals();
                    const auto &fun        = fe->get_phi();
                    const auto &JxW        = fe->get_JxW();

                    const uint n_fun = fun.size();
                    const uint n_qp  = fun[0].size();

                    vec.resize(n_fun * n_dims);
                    vec.zero();

                    for(uint qp = 0; qp < quad.n_points(); ++qp) {
                        for(uint i = 0; i < fun.size(); ++i){
                            for(uint d = 0; d < n_dims; ++d) {
                                vec(i + d * n_fun) += fun[i][qp] * fe_normals[qp](d) * JxW[qp];
                            }
                        }
                    }

                    dof_map.dof_indices(&e, dof_indices);

                    for(uint i = 0; i < dof_indices.size(); ++i) {
                        const uint ind = dof_indices[i];
                        for(uint d = 0; d < spatial_dim; ++d) {
                            normals.add(ind * spatial_dim + d, vec(i + d * n_fun));
                        }
                    }
                }
            }
        } //synch-block end

        auto r = range(normals);

        ReadAndWrite<UVector> rw(normals);
        for(auto i = r.begin(); i < r.end(); i += spatial_dim) {
            double len = 0.0;

            for(int d = 0; d < spatial_dim; ++d) {
                auto x = normals.get(i + d);
                len += x*x;
            }

            len = std::sqrt(len);

            if(len > 0.0) {
                for(int d = 0; d < spatial_dim; ++d) {
                    auto x = normals.get(i + d);
                    normals.set(i + d, x/len);
                }
            }
        }

        return true;
    }


    void CurvatureBasedEdgeProjection::read(Input &is)
    {
        std::cout << "CurvatureBasedEdgeProjection::read ..." << std::endl;

        is.get("sides", [this](Input &in) {
            in.get_all([this](Input &in) {
                int id = -1;
                in.get("id", id);
                boundary_tags.push_back(id);

                // std::cout << "side: " << id << std::endl;
            });
        });
    }

    static void make_index(const moonolith::Polygon<double, 3> &poly, const libMesh::Elem &elem, std::vector<uint> &idx)
    {
        auto n = poly.size();
        idx.resize(n);

        assert(n <= elem.n_nodes());

        moonolith::Vector<double, 3> p;
        for(std::size_t i = 0; i < n; ++i) {
            make(elem.node_ref(i), p);
            for(std::size_t j = 0; j < n; ++j) {
                if(distance(p, poly[j]) < 1e-16) {
                    idx[j] = i;
                    break;
                }
            }
        }
    }

    void CurvatureBasedEdgeProjection::apply(libMesh::MeshBase &mesh)
    {
        using Point2 = moonolith::Vector<double, 2>;
        using Point3 = moonolith::Vector<double, 3>;


        LibMeshFunctionSpace V(mesh, libMesh::LAGRANGE, libMesh::FIRST);
        V.initialize();

        UVector is_normal_component;
        UVector normals;
        USparseMatrix mat;
        // std::vector<int> boundary_tags;

        bool ok = assemble_smooth_normals(
            mesh,
            V.dof_map(),
            boundary_tags,
            normals
            ); assert(ok);

        UVector trafo_points = local_zeros(mesh.n_local_nodes()*3);
        UVector touched = local_zeros(mesh.n_local_nodes());

        auto vr = range(trafo_points);

        moonolith::PNTriangle<double, 3> pn_triangle;
        std::array<Point3, 3> triangle, t_normals;
        moonolith::Polygon<double, 3> polygon;
        std::vector<uint> node_idx;

        moonolith::AffineTransform<double, 2, 3> trafo;

        // int export_id = 0;
        // moonolith::MatlabScripter matlab;
        // matlab.close_all();



        {
            Read<UVector> r(normals);
            Write<UVector> w(trafo_points), wt(touched);
            for(auto e_it = mesh.active_local_elements_begin(); e_it != mesh.active_local_elements_end(); ++e_it) {
                const auto &e = **e_it;

                for(uint side = 0; side < e.n_sides(); ++side) {
                    if(e.neighbor_ptr(side) != nullptr) {continue;}

                    //if empty take them all
                    bool select = boundary_tags.empty();
                    for(auto t : boundary_tags) {
                        if(mesh.get_boundary_info().has_boundary_id(&e, side, t)) {
                            select = true;
                            break;
                        }
                    }

                    if(!select) continue;

                    auto side_ptr = e.build_side_ptr(side);
                    make(*side_ptr, polygon);
                    make_index(polygon, *side_ptr, node_idx);
                    make_transform(*side_ptr, trafo);

                    for(std::size_t i = 0; i < 3; ++i) {
                        const auto dof = side_ptr->node_ref(node_idx[i]).dof_number(0, 0, 0);
                        triangle[i] = polygon[i];

                        for(int d = 0; d < 3; ++d) {
                            t_normals[i][d] = normals.get(dof * 3 + d);
                        }
                    }

                    pn_triangle.init(triangle, t_normals);

                    // bool must_export = export_id == 0;

                    // if(must_export) {

                    // write(matlab,
                    //        triangle,
                    //        t_normals,
                    //        pn_triangle
                    // );
                    // }

                    // export_id++;

                    // exit(0);

                    for(int i = 0; i < side_ptr->n_nodes(); ++i) {
                        const auto dof = side_ptr->node_id(i);
                        if(!vr.inside(dof)) continue;

                        Point2 ref_point;
                        Point3 p;
                        make(side_ptr->node_ref(i), p);
                        trafo.apply_inverse(p, ref_point);


                        pn_triangle.apply(ref_point, p);

                        // if(must_export) {
                        //     print(ref_point);
                        //     matlab.plot(p, "b.");
                        //     matlab.text(p, std::to_string(i));
                        // }

                        touched.set(dof, 1.);
                        for(int d = 0; d < 3; ++d) {
                            trafo_points.set(dof * 3 + d, p[d]);

                        }
                    }
                }
            }
        }


        {
            Read<UVector> r(trafo_points), rt(touched);
            // UVector g_points = ghosted

            for(auto n_it = mesh.local_nodes_begin(); n_it != mesh.local_nodes_end(); ++n_it) {
                auto &node = **n_it;
                const auto dof = node.id();

                if(touched.get(dof) > 0.) {
                    for(int d = 0; d < 3; ++d) {
                       node(d) = trafo_points.get(dof * 3 + d);
                    }
                }
            }
        }

        // matlab.axis_equal();
        // matlab.save("mesh.m");
    }
}
