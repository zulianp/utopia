//https://fossies.org/diffs/petsc/3.7.7_vs_3.8.0/
//https://www.mcs.anl.gov/petsc/documentation/changes/39.html

#include "utopia_petsc.hpp"
#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_RMTR.hpp"
#include "utopia_RMTR_inf.hpp"
#include "utopia_MG_OPT.hpp"
#include "utopia_FAS.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_AffineSimilarity.hpp"

//explicit instantiations
namespace utopia {
 	template class Wrapper<PetscSparseMatrix, 2>;
	template class Wrapper<PetscMatrix, 2>;
	template class Wrapper<PetscVector, 1>;

	//petsc linear solvers and smoothers
	template class KSPSolver<DSMatrixd, DVectord>;
	template class ConjugateGradient<DSMatrixd, DVectord>;
	template class GaussSeidel<DSMatrixd, DVectord>;
	

	//petsc non-linear solvers
	template class NonLinearGaussSeidel<DSMatrixd, DVectord>;
	
	template class Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL>;
	template class RMTR<DSMatrixd, DVectord, FIRST_ORDER>;
	template class RMTR_inf<DSMatrixd, DVectord>;

	template class FAS<DSMatrixd, DVectord>;
	template class MG_OPT<DSMatrixd, DVectord>;

	template class AffineSimilarity<DSMatrixd, DVectord>;

	void optimize_nnz(DSMatrixd &A)
	{
		auto rr = row_range(A);
		auto cr = A.implementation().col_range();
		auto ls = local_size(A);
		auto gs = size(A);

		std::vector<PetscInt> d_nnz(rr.extent(), 0), o_nnz(rr.extent(), 0);
		each_read(A, [&](const utopia::SizeType i, const utopia::SizeType j, const PetscScalar val) {
			if(std::abs(val) > 1e-18) {
				if(cr.inside(j)) {
					++d_nnz[i - rr.begin()];
				} else {
					++o_nnz[i - rr.begin()];
				}
			}
		});

		DSMatrixd A_opt;

		A_opt.implementation().matij_init(
			A.implementation().communicator(),
			A.implementation().type_override(),
			ls.get(0),
			ls.get(1),
			gs.get(0),
			gs.get(1),
			d_nnz,
			o_nnz
			);

		{
			Write<DSMatrixd> w_A(A_opt);
			each_read(A, [&](const SizeType i, const SizeType j, const PetscScalar val) {
				if(std::abs(val) > 1e-18) {
					A_opt.set(i, j, val);
				}
			});

		}

		A = std::move(A_opt);
	}

	bool is_diagonally_dominant(const DSMatrixd &A)
	{
		DVectord d = diag(A);
		DVectord o = local_zeros(local_size(d));

		{
			Write<DVectord> w_o(o);
			each_read(A,[&o](const SizeType i, const SizeType j, const PetscScalar val) {
				if(i != j) {
					o.add(i, std::abs(val));
				}
			});
		}

		DVectord diff = d - o;
		PetscScalar m = min(diff);
		return m > 0.;
	}
}

