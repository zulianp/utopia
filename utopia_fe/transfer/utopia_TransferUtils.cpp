
#include "utopia_TransferUtils.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"

#include <memory>

namespace utopia {

	bool assemble_interpolation(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D)
	{
	    auto assembler = std::make_shared<InterpolationLocalAssembler>(from.mesh().mesh_dimension());
	    auto local2global = std::make_shared<Local2Global>(true);

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
	    D = diag(sum(B, 1));

	    double sum_B = sum(B);
	    double sum_D = sum(D);

	    std::cout << sum_B << " == " << sum_D << std::endl;
	    return true;
	}

	bool assemble_projection(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D, const bool use_biorth)
	{
	    bool is_shell = from.mesh().mesh_dimension() < from.mesh().spatial_dimension();

	    auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), use_biorth, true, is_shell);
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

}
