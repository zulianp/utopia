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
        std::vector<int> boundary_tags;

        bool ok = assemble_smooth_normals(
            mesh,
            V.dof_map(),
            boundary_tags,
            normals
        ); assert(ok);




        moonolith::PNTriangle<double, 3> pn_triangle;
        std::array<Point3, 3> triangle, edge_points, t_normals;
        std::array<Point2, 3> edge_nodes
        { 
            Point2(0.5, 0.0),
            Point2(0.5, 0.5),
            Point2(0.0, 0.5)
        };


        Read<UVector> r(normals);
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

                auto side_ptr = e.build_side_ptr(side);

                int vert = 0;
                for(int i = 0; i < side_ptr->n_nodes(); ++i) {
                    if(side_ptr->is_vertex(i)) {
                        make(side_ptr->node_ref(i), triangle[vert]);
                        const auto dof = side_ptr->node_ref(i).dof_number(0, 0, 0);

                        for(int d = 0; d < 3; ++d) {
                            t_normals[vert][d] = normals.get(dof * 3 + d);
                        }

                        vert++;
                    }
                }

                pn_triangle.init(triangle, t_normals);

                std::cout << "===========================\n";

                for(int j = 0; j < 3; ++j) {
                    print(triangle[j]);
                }

                std::cout << "---------------------------\n";
                for(int j = 0; j < 3; ++j) {

                    pn_triangle.apply(edge_nodes[j], edge_points[j]);
                    print(edge_points[j]);
                }

                std::cout << "===========================\n";


                vert = 0;
                for(int i = 0; i < side_ptr->n_nodes(); ++i) {
                    if(!side_ptr->is_vertex(i)) {
                        auto &node = mesh.node_ref(side_ptr->node(i));

                        for(int d = 0; d < 3; ++d) {
                            node(d) = edge_points[vert][d];
                        }

                        ++vert;
                    }
                }
            }
        }
    }
}
