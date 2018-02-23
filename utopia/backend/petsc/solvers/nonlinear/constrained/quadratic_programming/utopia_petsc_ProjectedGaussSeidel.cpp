#include "utopia_petsc_ProjectedGaussSeidel.hpp"
#include "utopia_petsc_external.hpp"

extern PetscErrorCode constrain_box(MatScalar * diag, PetscScalar * t, PetscScalar * x,
                              PetscScalar * lb, PetscScalar * ub, int row,
                              int bs, PetscBool * constrained);

extern PetscErrorCode constrain_blocknd(MatScalar * diag, PetscScalar * t, PetscScalar * x,
                                  PetscScalar * lb, PetscScalar * ub, int row,
                                  int bs, PetscBool * constrained);

namespace utopia {
	
	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::solve(const Matrix &A, const Vector &b, Vector &x)
	{
		init(A);

		auto ctx = std::static_pointer_cast<NBGS_CTX>(ctx_);

		if(empty(x) || size(x).get(0) != size(b).get(0)) {
			x = local_zeros(local_size(b));
		}

		old_x_ = x;

		Vector ub, lb;

		if(constraints_.has_lower_bound()) {
			lb = *constraints_.lower_bound();
		} else {
			lb = local_values(local_size(b).get(0), -std::numeric_limits<PetscScalar>::max());
		}

		if(constraints_.has_upper_bound()) {
			ub = *constraints_.upper_bound();
		} else {
			ub = local_values(local_size(b).get(0), std::numeric_limits<PetscScalar>::max());
		}

		ctx->ub = raw_type(ub);
		ctx->lb = raw_type(lb);

		for(SizeType i = 0; i < this->max_it(); ++i) {
			NBGS(ctx.get(), raw_type(A), raw_type(b), raw_type(x), raw_type(lb), raw_type(ub), ctx->blocksize);

			PetscScalar diff = norm2(old_x_ - x);
			if(approxeq(diff, 0., this->stol())) {
				std::cout << "CONVERGED" << std::endl;
				return true;
			}

			prev_active_set_ = active_set_;

			if(this->verbose()) {
				std::cout << "iter: " << i << " diff: " << diff << std::endl;
			}

			old_x_ = x;
		}

		return false;
	}

	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::smooth(const Matrix &A, const Vector &b, Vector &x)
	{
		std::cerr << "[Error] not implemented" << std::endl;
		return false;
	}

	ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::~ProjectedGaussSeidel() {}

	void ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::init(const Matrix &A)
	{
		auto ls = local_size(A);
		auto gs = size(A);

		active_set_      = local_zeros(ls.get(0));
		prev_active_set_ = local_zeros(ls.get(0));
		changed_         = local_zeros(ls.get(0));

		auto ctx = std::make_shared<NBGS_CTX>();

		ctx->active_coordinates      = raw_type(active_set_);
		ctx->prev_active_coordinates = raw_type(prev_active_set_);
		ctx->changed_coordinates     = raw_type(changed_);
		
		ctx->alpha = 1.;
		ctx->blocksize = block_size_;
		ctx->constrained = constraints_.has_bound()? PETSC_TRUE : PETSC_FALSE;
		ctx->with_linesearch = use_line_search_? PETSC_TRUE : PETSC_FALSE;


		if(block_size_ == 1) {
			ctx->constrain = &constrain_box;
		} else {
			ctx->constrain = &constrain_blocknd;
		}

		ctx_ = ctx;
	}
}
