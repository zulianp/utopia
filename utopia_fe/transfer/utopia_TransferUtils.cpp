
#include "utopia_TransferUtils.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_QuadratureUtils.hpp"


#include <memory>
#include <cassert>
#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/dof_map.h>

namespace utopia {

    bool assemble_interpolation(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D, const int n_var)
    {
        auto assembler = std::make_shared<InterpolationLocalAssembler>(from.mesh().mesh_dimension());
        auto local2global = std::make_shared<Local2Global>(true);


        TransferOptions opts;
        opts.n_var        = n_var;

        TransferAssembler transfer_assembler(assembler, local2global);

        std::vector< std::shared_ptr<USparseMatrix> > mats;
        if(!transfer_assembler.assemble(
                                        make_ref(from.mesh()),
                                        make_ref(from.dof_map()),
                                        make_ref(to.mesh()),
                                        make_ref(to.dof_map()),
                                        mats,
                                        opts)) {
            return false;
        }

        B = std::move(*mats[0]);
        D = diag(sum(B, 1));

        double sum_B = sum(B);
        double sum_D = sum(D);

        std::cout << sum_B << " == " << sum_D << std::endl;
        return true;
    }

    bool assemble_projection(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D, const bool use_biorth, const int n_var)
    {
        bool is_shell = from.mesh().mesh_dimension() < from.mesh().spatial_dimension();

        auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), use_biorth, true, is_shell);
        auto local2global = std::make_shared<Local2Global>(false);

        TransferAssembler transfer_assembler(assembler, local2global);

        TransferOptions opts;
        opts.n_var        = n_var;

        std::vector< std::shared_ptr<USparseMatrix> > mats;
        if(!transfer_assembler.assemble(
                                        make_ref(from.mesh()),
                                        make_ref(from.dof_map()),
                                        make_ref(to.mesh()),
                                        make_ref(to.dof_map()),
                                        mats,
                                        opts)) {
            return false;
        }


        B = std::move(*mats[0]);

        if(use_biorth) {
            UVector d = sum(B, 1);
      //   	UVector d_inv = local_values(local_size(d), 1.);

            // {
            // 	Write<UVector> w(d_inv);
         //    	each_read(d, [&d_inv](const SizeType i, const double value) {
         //    		if(std::abs(value) > 1e-14) {
         //    			d_inv.set(i, 1./value);
         //    		}
         //    	});
      //   	}

            D = diag(d);

        } else {

            D = std::move(*mats[1]);
        }

        double sum_B = sum(B);
        double sum_D = sum(D);

        std::cout << sum_B << " == " << sum_D << std::endl;
        return true;
    }

    bool assemble_coupling(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B)
    {
        bool is_shell = from.mesh().mesh_dimension() < from.mesh().spatial_dimension();

        auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), false, false, is_shell);
        auto local2global = std::make_shared<Local2Global>(false);

        TransferAssembler transfer_assembler(assembler, local2global);

        std::vector< std::shared_ptr<USparseMatrix> > mats;
        if(!transfer_assembler.assemble(
                                        make_ref(from.mesh()),
                                        make_ref(from.dof_map()),
                                        make_ref(to.mesh()),
                                        make_ref(to.dof_map()),
                                        mats)) {
            return false;
        }

        B = std::move(*mats[0]);

        double sum_B = sum(B);

        std::cout << sum_B << std::endl;
        return true;
    }

    //use different lagr mult space
    bool assemble_projection(
        LibMeshFunctionSpace &from,
        LibMeshFunctionSpace &to,
        LibMeshFunctionSpace &lagr,
        USparseMatrix &B, USparseMatrix &D)
    {
        if(assemble_coupling(from, lagr, B)) {
            return assemble_coupling(to, lagr, D);
        } else {
            return false;
        }
    }

    bool assemble_interpolation(
        const libMesh::MeshBase &mesh, 
        const libMesh::DofMap &dof_map_p1,
        const libMesh::DofMap &dof_map_p2,
        USparseMatrix &mat
    )
    {
        auto nnz = max_nnz_x_row(dof_map_p2);

        mat = local_sparse(dof_map_p2.n_local_dofs(), dof_map_p1.n_local_dofs(), nnz);
        Write<USparseMatrix> w(mat, utopia::GLOBAL_INSERT);

        const auto dim = mesh.mesh_dimension();
       
        auto sys_p1 = dof_map_p1.sys_number();
        auto sys_p2 = dof_map_p2.sys_number();
        
        auto n_vars = dof_map_p2.n_variables();

        assert( n_vars == dof_map_p1.n_variables() );

        auto fe_p1 = libMesh::FEBase::build(dim, dof_map_p1.variable_type(0));
        auto fe_p2 = libMesh::FEBase::build(dim, dof_map_p2.variable_type(0));

        const auto &phi1 = fe_p1->get_phi();
        const auto &phi2 = fe_p2->get_phi();

        libMesh::DenseMatrix<double> el_mat;
        std::vector<libMesh::dof_id_type> dof_p1, dof_p2;

        for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
            const auto &e = **e_it;

            auto q = QuadratureUtils::nodal_quad_points(e);

            fe_p1->attach_quadrature_rule(q.get());
            fe_p2->attach_quadrature_rule(q.get());

            fe_p1->reinit(&e);
            fe_p2->reinit(&e);

            for(uint var = 0; var < n_vars; ++var) {
                dof_map_p1.dof_indices(&e, dof_p1, var);
                dof_map_p2.dof_indices(&e, dof_p2, var);

                const auto n_dofs_p1 = dof_p1.size();
                const auto n_dofs_p2 = dof_p2.size();

                el_mat.resize(n_dofs_p2, n_dofs_p1);
                el_mat.zero();

                const auto n_qp = phi1[0].size();
                const auto n_i  = phi2.size();
                const auto n_j  = phi1.size();

                for(uint i = 0; i < n_i; ++i) {
                    for(uint j = 0; j < n_j; ++j) {
                        for(uint k = 0; k < n_qp; ++k) {
                            auto val = phi2[i][k] * phi1[j][k];
                            if(std::abs(val) < 1e-10) continue;
                            
                            el_mat(i, j) = val;
                        }
                    }
                }

                //local 2 global
                for(uint i = 0; i < n_dofs_p2; ++i) {
                    for(uint j = 0; j < n_dofs_p1; ++j) {
                        auto val = el_mat(i, j);
                        if(std::abs(val) < 1e-10) continue;

                        mat.set_matrix({dof_p2[i]}, {dof_p1[j]}, {val});
                    }
                }
            }
        }

        return true;
    }

}
