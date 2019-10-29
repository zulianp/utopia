#ifndef UTOPIA_LIBMESH_ASSEMBLE_LOCAL_HPP
#define UTOPIA_LIBMESH_ASSEMBLE_LOCAL_HPP

#include "utopia_FEEval_Local.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FiniteElement.hpp"


namespace utopia {



    template<class Derived, class FE, class LocalAssembler, class InitExpr>
    void libmesh_assemble_aux(
        FunctionSpace<Derived> &V_w,
        FE &element,
        const Expression<InitExpr> &init_expr,
        USparseMatrix &matrix,
        LocalAssembler assembler
        )
    {
        const bool disable_adaptivity = utopia::Utopia::instance().get("disable-adaptivity") == "true";

        Chrono c;
        c.start();

        auto &V = V_w.derived();
        auto &mesh = V.mesh();
        auto &dof_map = V.dof_map();

       
        libMesh::DofConstraints constraints;

        if(!disable_adaptivity) {
            Adaptivity::compute_all_constraints(
                mesh,
                dof_map,
                constraints
            );
        }

        LibMeshAssembler::allocate_matrix(
            V.dof_map(),
            matrix
        );


        USerialMatrix mat;

        std::vector<libMesh::dof_id_type> dofs;

        {
            Write<USparseMatrix> w(matrix, utopia::GLOBAL_ADD);

            auto e_begin = mesh.active_local_elements_begin();

            //Which basis-function
            auto u = trial(element);
            for(auto e_it = e_begin; e_it != mesh.active_local_elements_end(); ++e_it) {
                element.set((*e_it)->id());

                if(e_it == e_begin) {
                    element.init(init_expr.derived());
                } else {
                    element.reinit(init_expr.derived());
                }

                assembler(element, mat);

                // if(element.ctx().has_assembled()) {
                    V.dof_map().dof_indices(*e_it, dofs);

                    if(!disable_adaptivity) {
                        Adaptivity::constrain_matrix(*e_it, dof_map, constraints, mat, dofs);
                    }


                    add_matrix(mat, dofs, dofs, matrix);
                // }
            }
        }

        c.stop();
        std::cout << "assemble mat(E=" << mesh.n_active_local_elem() << "): " << c << std::endl;
    }

    template<class Derived, class FE, class LocalAssembler, class InitExpr>
    void libmesh_assemble_aux(
        FunctionSpace<Derived> &V_w,
        FE &element,
        const Expression<InitExpr> &init_expr,
        UVector &vector,
        LocalAssembler assembler)
    {
        const bool disable_adaptivity = utopia::Utopia::instance().get("disable-adaptivity") == "true";

        Chrono c;
        c.start();

        auto &V = V_w.derived();
        auto &mesh = V.mesh();
        auto &dof_map = V.dof_map();

        libMesh::DofConstraints constraints;

        if(!disable_adaptivity) {
            Adaptivity::compute_all_constraints(
                mesh,
                dof_map,
                constraints
            );
        }

        typename Traits<UVector>::IndexSet ghost_nodes;
        convert(dof_map.get_send_list(), ghost_nodes);

        UVector temp_vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);

        {
            Write<UVector> w(temp_vec, utopia::GLOBAL_ADD);

            USerialVector vec;

            std::vector<libMesh::dof_id_type> dofs;
        

            auto e_begin = mesh.active_local_elements_begin();

            //Which basis-function
            auto u = trial(element);
            for(auto e_it = e_begin; e_it != mesh.active_local_elements_end(); ++e_it) {
                element.set((*e_it)->id());

                if(e_it == e_begin) {
                    element.init(init_expr.derived());
                } else {
                    element.reinit(init_expr.derived());
                }

                assembler(element, vec);

                // if(element.ctx().has_assembled()) {
                    V.dof_map().dof_indices(*e_it, dofs);

                    if(!disable_adaptivity) {
                        Adaptivity::constrain_vector(*e_it, dof_map, constraints, vec, dofs);
                    }

                    add_vector(vec, dofs, temp_vec);
                // }
            }
        }

        if(Traits<UVector>::Backend == utopia::TRILINOS) {
            vector = 1. * temp_vec;
        } else {
            vector = std::move(temp_vec);
        }

        c.stop();
        std::cout << "assemble vec(E=" << mesh.n_active_local_elem() << "): " << c << std::endl;
    }


    template<class LocalAssembler, class InitExpr>
    void assemble(
        FunctionSpace<LibMeshFunctionSpace> &V_w,
        const Expression<InitExpr> &init_expr,
        USparseMatrix &matrix,
        LocalAssembler assembler
        )
    {
        FiniteElement<LibMeshFunctionSpace> element(V_w.derived());
        libmesh_assemble_aux(V_w, element, init_expr, matrix, assembler);
    }

    template<class LocalAssembler, class InitExpr>
    void assemble(
        ProductFunctionSpace<LibMeshFunctionSpace> &V_w,
        const Expression<InitExpr> &init_expr,
        USparseMatrix &matrix,
        LocalAssembler assembler
        )
    {
        FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace> > element(V_w.derived());
        libmesh_assemble_aux(V_w.derived().subspace(0), element, init_expr, matrix, assembler);
    }


    template<class LocalAssembler, class InitExpr>
    void assemble(
        FunctionSpace<LibMeshFunctionSpace> &V_w,
        const Expression<InitExpr> &init_expr,
        UVector &vector,
        LocalAssembler assembler
        )
    {
        FiniteElement<LibMeshFunctionSpace> element(V_w.derived());
        libmesh_assemble_aux(V_w, element, init_expr, vector, assembler);
    }

    template<class LocalAssembler, class InitExpr>
    void assemble(
        ProductFunctionSpace<LibMeshFunctionSpace> &V_w,
        const Expression<InitExpr> &init_expr,
        UVector &vector,
        LocalAssembler assembler
        )
    {
        FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace> > element(V_w.derived());
        libmesh_assemble_aux(V_w.derived().subspace(0), element, init_expr, vector, assembler);
    }

}

#endif //UTOPIA_LIBMESH_ASSEMBLE_LOCAL_HPP