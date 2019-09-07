#ifndef UTOPIA_UI_MESH_HPP
#define UTOPIA_UI_MESH_HPP

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_modification.h"
#include "utopia_UIMorph.hpp"
#include "utopia_MeshParamSmoother.hpp"

#include <memory>

namespace utopia {

    inline void refine_0(libMesh::UnstructuredMesh &mesh)
    {
        libMesh::MeshRefinement mesh_refinement(mesh);
        auto e_it = elements_begin(mesh);

        if(e_it != elements_end(mesh)) {
            (*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
        }

        mesh_refinement.make_flags_parallel_consistent();
        mesh_refinement.refine_elements();
        mesh_refinement.test_level_one(true);

        mesh.prepare_for_use();
    }

    inline void random_refine(
        libMesh::UnstructuredMesh &mesh,
        const int refinement_loops
        )
    {

        libMesh::MeshRefinement mesh_refinement(mesh);

        for(int i = 0; i < refinement_loops; ++i) {
            mesh_refinement.clean_refinement_flags();

            int idx = 0;
            for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it, ++idx) {
                auto val = idx % 2 == 1;
                if(val) {
                    (*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
                }
            }

            mesh_refinement.make_flags_parallel_consistent();
            mesh_refinement.refine_elements();
            mesh_refinement.test_level_one(true);
        }

        // mesh_refinement.clean_refinement_flags();
        mesh.prepare_for_use();
    }

    template<class Mesh>
    class UIMesh {};// final : public Configurable { };

    template<>
    class UIMesh<libMesh::DistributedMesh> : public Configurable {
    public:
        template<class... Args>
        UIMesh(Args &&...args)
        : mesh_(std::make_shared<libMesh::DistributedMesh>(std::forward<Args...>(args...))), empty_(true)
        {}

        void read(Input &is) override {
            std::string mesh_type = "square";
            std::string path = "";

            empty_ = false;

            int refinements = 0;

            double span[3] = { 0., 0., 0. };

            double min_coords[3] = {0., 0., 0.};
            double max_coords[3] = {1., 1., 1.};

            int n[3] = {5, 5, 5};

            double scale = 1.;
            double shift[3] = {0. , 0., 0.};

            bool full_order = true;

            std::string elem_type = "quad";

            is.get("type", mesh_type);
            is.get("elem-type", elem_type);
            is.get("order", order);
            is.get("full-order", full_order);
            is.get("path", path);

            is.get("refinements", refinements);

            is.get("span-x", span[0]);
            is.get("span-y", span[1]);
            is.get("span-z", span[2]);

            is.get("min-x", min_coords[0]);
            is.get("min-y", min_coords[1]);
            is.get("min-z", min_coords[2]);

            is.get("max-x", max_coords[0]);
            is.get("max-y", max_coords[1]);
            is.get("max-z", max_coords[2]);

            is.get("n-x", n[0]);
            is.get("n-y", n[1]);
            is.get("n-z", n[2]);

            is.get("scale", scale);
            is.get("shift-x", shift[0]);
            is.get("shift-y", shift[1]);
            is.get("shift-z", shift[2]);


            if(mesh_type == "file") {
                mesh_->read(path);
            } else if(mesh_type == "line") {
                libMesh::MeshTools::Generation::build_line(
                    *mesh_, n[0], min_coords[0], max_coords[0], get_type(elem_type, order, 1)
                );
            } else if(mesh_type == "square") {
                libMesh::MeshTools::Generation::build_square(*mesh_,
                    n[0], n[1],
                    min_coords[0], max_coords[0],
                    min_coords[1], max_coords[1],
                    get_type(elem_type, order, 2)
                    );
            } else if(mesh_type == "cube") {
                libMesh::MeshTools::Generation::build_cube(*mesh_,
                    n[0], n[1], n[2],
                    min_coords[0], max_coords[0],
                    min_coords[1], max_coords[1],
                    min_coords[2], max_coords[2],
                    get_type(elem_type, order, 3)
                    );
            } else if(mesh_type == "sphere") {

                double radius = 1.;
                int sphere_refine = 2;

                is.get("radius", radius);
                is.get("sphere-refine", sphere_refine);

                libMesh::MeshTools::Generation::build_sphere(*mesh_,
                    radius,
                    sphere_refine,//const unsigned int nr = 2,
                    get_type(elem_type, order, 3)
                    //const unsigned int 	n_smooth = 2,
                    // const bool 	flat = true
                );

            } else if(mesh_type == "aabb") {
                libMesh::DistributedMesh temp_mesh(mesh_->comm());
                temp_mesh.read(path);

                auto bb = bounding_box(temp_mesh);

                if(temp_mesh.spatial_dimension() == 3) {
                    libMesh::MeshTools::Generation::build_cube(
                        *mesh_,
                        n[0], n[1], n[2],
                        bb.min()(0) - span[0], bb.max()(0) + span[0],
                        bb.min()(1) - span[1], bb.max()(1) + span[1],
                        bb.min()(2) - span[2], bb.max()(2) + span[2],
                        get_type(elem_type, order, 3, full_order)
                        );

                } else {
                    libMesh::MeshTools::Generation::build_square(
                        *mesh_,
                        n[0], n[1],
                        bb.min()(0) - span[0], bb.max()(0) + span[0],
                        bb.min()(1) - span[1], bb.max()(1) + span[1],
                        get_type(elem_type, order, 2, full_order)
                        );
                }
            }

            int block_override = -1;
            is.get("block-override", block_override);

            override_block(block_override, *mesh_);

            //build_extrusion (UnstructuredMesh &mesh, const MeshBase &cross_section, const unsigned int nz, RealVectorValue extrusion_vector, QueryElemSubdomainIDBase *elem_subdomain=libmesh_nullptr)

            scale_mesh(scale, *mesh_);
            shift_mesh(shift, *mesh_);

            refine(refinements, *mesh_);

            bool convert_to_triangles = false;
            is.get("convert-to-triangles", convert_to_triangles);

            if(convert_to_triangles) {
                libMesh::MeshTools::Modification::all_tri(*mesh_);
            }

            if(mesh_type == "file" && order == 2) {
                mesh_->all_second_order(full_order);
            }

            { 
                //single morph
                UIMorph<libMesh::DistributedMesh> morph;
                is.get("morph", morph);
                
                if(morph.is_valid()) {
                    morph.apply(*mesh_);
                }
            }

            // is.get("improve-smoothness", [this](Input &is) {
            //     MeshParamSmoother smoother;
            //     smoother.read(is);
            //     smoother.apply(*mesh_);
            // });

            //multi morph
            is.get("morphs", [this](Input &is) {
                is.get_all([this](Input &is) {

                    UIMorph<libMesh::DistributedMesh> morph;
                    morph.read(is);
                    
                    if(morph.is_valid()) {
                        morph.apply(*mesh_);
                    }
                    
                });
            });


            bool must_refine_0 = false;
            is.get("refine-0", must_refine_0);

            bool must_random_refine = false;
            is.get("random-refine", must_random_refine);

            int n_random_refinements = 1;
            is.get("random-refinements", n_random_refinements);

            if(must_refine_0) {
                refine_0(*mesh_);
            }

            if(must_random_refine) {
                random_refine(*mesh_, n_random_refinements);
            }

        }

        inline libMesh::DistributedMesh &mesh()
        {
            assert(mesh_);
            return *mesh_;
        }

        inline std::shared_ptr<libMesh::DistributedMesh> mesh_ptr()
        {
            return mesh_;
        }

        inline bool empty() const {
            return empty_;
        }

    private:
        int order = 1;
        std::shared_ptr<libMesh::DistributedMesh> mesh_;
        bool empty_;

        ////////////////////////////////////////////////

        libMesh::ElemType get_type(
            const std::string &elem_type,
            const int order,
            const int dim,
            const bool full_order = true) const
        {
            if(dim == 3) {
                libMesh::ElemType type = libMesh::HEX8;

                if(elem_type == "tet") {
                    type = libMesh::TET4;
                }

                if(order == 2) {
                    type = libMesh::HEX20;

                    if(elem_type == "tet") {
                        type = libMesh::TET10;
                    } else if(full_order) {
                        type = libMesh::HEX27;
                    }
                }

                if(order == 3) {
                    type = libMesh::HEX27;
                }

                if(elem_type == "prism") {
                    type = libMesh::PRISM6;

                    if(order == 2) {
                        type = libMesh::PRISM15;
                    }
                }

                if(elem_type == "pyramid") {
                    type = libMesh::PYRAMID5;

                    if(order == 2) {
                        type = libMesh::PYRAMID13;
                    }
                }

                return type;

            } else if(dim == 2) {
                libMesh::ElemType type = libMesh::QUAD4;

                if(elem_type == "tri") {
                    type = libMesh::TRI3;
                }

                if(order == 2) {
                    type = libMesh::QUAD8;

                    if(elem_type == "tri") {
                        type = libMesh::TRI6;
                    }
                }

                if(order == 3) {
                    type = libMesh::QUAD9;

                     if(elem_type == "tri") {
                        assert(false);
                        // type = libMesh::TRI6;
                    }
                }

                return type;
            } else if(dim == 1) {
                libMesh::ElemType type = libMesh::EDGE2;

                if(order == 2) {
                    type = libMesh::EDGE3;
                }

                return type;
            }

            return libMesh::TRI3;
        }

        static void refine(const int n_refs, libMesh::MeshBase &mesh)
        {
            if(n_refs <= 0) return;

            libMesh::MeshRefinement mesh_refinement(mesh);
            mesh_refinement.make_flags_parallel_consistent();
            mesh_refinement.uniformly_refine(n_refs);
        }

        static void scale_mesh(const double &scale_factor, libMesh::MeshBase &mesh)
        {
            if(scale_factor == 1.) return;
            assert(scale_factor > 0.);

            for(auto it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it) {

                for(int i = 0; i < LIBMESH_DIM; ++i) {
                    (**it)(i) *= scale_factor;
                }
            }
        }


        static void override_block(const int block, libMesh::MeshBase &mesh)
        {
            if(block <= 0) return;

            for(auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
                (*it)->subdomain_id() = block;
            }
        }


        static void shift_mesh(const double t[3], libMesh::MeshBase &mesh)
        {
            if(0. == t[0] && 0. == t[1] && 0. == t[2]) return;

            for(auto it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it) {
                for(int i = 0; i < LIBMESH_DIM; ++i) {
                    (**it)(i) += t[i];
                }
            }
        }
    };
}


#endif //UTOPIA_UI_MESH_HPP
