#include "utopia_petsc_ProjectedGaussSeidel.hpp"
#include "utopia_petsc_external.hpp"

extern PetscErrorCode constrain_box(MatScalar * diag, PetscScalar * t, PetscScalar * x,
	PetscScalar * lb, PetscScalar * ub, int row,
	int bs, PetscBool * constrained);

extern PetscErrorCode constrain_blocknd(MatScalar * diag, PetscScalar * t, PetscScalar * x,
	PetscScalar * lb, PetscScalar * ub, int row,
	int bs, PetscBool * constrained);

namespace utopia {
	
	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::apply(const Vector &b, Vector &x)
	{
		if(this->verbose())
			this->init_solver("utopia ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});

		const auto &A = *get_operator();

		auto ctx = std::static_pointer_cast<NBGS_CTX>(ctx_);

		if(empty(x) || size(x).get(0) != size(b).get(0)) {
			x = local_zeros(local_size(b));
		}

		old_x_ = x;
		const Vector &ub = *constraints_.upper_bound();
		const Vector &lb = *constraints_.lower_bound();

		ctx->ub = raw_type(ub);
		ctx->lb = raw_type(lb);

		bool converged = false;
		for(SizeType i = 0; i < this->max_it(); ++i) {
			NBGS(ctx.get(), raw_type(A), raw_type(b), raw_type(x), raw_type(lb), raw_type(ub), ctx->blocksize);

			PetscScalar diff = norm2(old_x_ - x);
			if(this->verbose())
				PrintInfo::print_iter_status({static_cast<Scalar>(i), diff}); 

			if(this->check_convergence(i, 1, 1, diff)) {
				converged = true;
				break;
			}

			prev_active_set_ = active_set_;
			old_x_ = x;
		}

		return converged;
	}

	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::smooth(const Matrix &A, const Vector &b, Vector &x)
	{
		if(empty(x) || size(x).get(0) != size(b).get(0)) {
			x = local_zeros(local_size(b));
		}

		init(A);

		const Vector &ub = *constraints_.upper_bound();
		const Vector &lb = *constraints_.lower_bound();

		auto ctx = std::static_pointer_cast<NBGS_CTX>(ctx_);
		ctx->ub = raw_type(ub);
		ctx->lb = raw_type(lb);

		for(SizeType i = 0; i < this->sweeps(); ++i) {
			NBGS(ctx.get(), raw_type(A), raw_type(b), raw_type(x), raw_type(lb), raw_type(ub), ctx->blocksize);
		}

		return true;
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


		if(!constraints_.has_lower_bound()) {
			constraints_.lower_bound() = std::make_shared<Vector>(local_values(ls.get(0), -std::numeric_limits<PetscScalar>::max()));
		}

		if(!constraints_.has_upper_bound()) {
			constraints_.upper_bound() = std::make_shared<Vector>(local_values(ls.get(0), std::numeric_limits<PetscScalar>::max()));
		}
	}
}
