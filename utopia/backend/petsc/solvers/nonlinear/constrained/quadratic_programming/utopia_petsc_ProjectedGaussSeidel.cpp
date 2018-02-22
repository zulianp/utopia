// #include "utopia_petsc_ProjectedGaussSeidel.hpp"
// #include "utopia_petsc_external.hpp"

// extern PetscErrorCode constrain_box(MatScalar * diag, PetscScalar * t, PetscScalar * x,
//                               PetscScalar * lb, PetscScalar * ub, int row,
//                               int bs, PetscBool * constrained);

// extern PetscErrorCode constrain_blocknd(MatScalar * diag, PetscScalar * t, PetscScalar * x,
//                                   PetscScalar * lb, PetscScalar * ub, int row,
//                                   int bs, PetscBool * constrained);

// namespace utopia {
	
// 	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::solve(const Matrix &A, const Vector &b, Vector &x)
// 	{
// 		init(A);

// 		auto ctx = std::static_pointer_cast<NBGS_CTX>(ctx_);

// 		if(empty(x) || size(x).get(0) != size(b).get(0)) {
// 			x = local_zeros(local_size(b));
// 		}

// 		Vector ub, lb;

// 		if(constraints_.has_lower_bound()) {
// 			lb = *constraints_.lower_bound();
// 		} else {
// 			lb = local_values(local_size(b).get(0), -std::numeric_limits<PetscScalar>::max());
// 		}

// 		if(constraints_.has_upper_bound()) {
// 			ub = *constraints_.upper_bound();
// 		} else {
// 			ub = local_values(local_size(b).get(0), std::numeric_limits<PetscScalar>::max());
// 		}

// 		ctx->ub = raw_type(ub);
// 		ctx->lb = raw_type(lb);

// 		for(SizeType i = 0; i < this->max_it(); ++i) {
// 			NBGS(ctx.get(), raw_type(A), raw_type(b), raw_type(x), raw_type(lb), raw_type(ub), ctx->blocksize);
// 			changed_ = abs(prev_active_set_ - active_set_);
// 			SizeType n_changed = double(sum(changed_));
// 			if(n_changed == 0) {
// 				return true;
// 			}

// 			prev_active_set_ = active_set_;

// 			if(this->verbose()) {
// 				std::cout << "iter: " << i << std::endl;
// 			}
// 		}

// 		return false;
// 	}

// 	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::smooth(const Matrix &A, const Vector &b, Vector &x)
// 	{
// 		std::cerr << "[Error] not implemented" << std::endl;
// 		return false;
// 	}

// 	ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::~ProjectedGaussSeidel() {}

// 	void ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::init(const Matrix &A)
// 	{
// 		auto ls = local_size(A);
// 		auto gs = size(A);

// 		// d_inv_ = diag(1./diag(A));

// 		active_set_      = local_zeros(ls.get(0));
// 		prev_active_set_ = local_zeros(ls.get(0));
// 		changed_         = local_zeros(ls.get(0));

// 		auto ctx = std::make_shared<NBGS_CTX>();
// 		// ctx->Dinv = &raw_type(d_inv_);

// 		ctx->active_coordinates      = raw_type(active_set_);
// 		ctx->prev_active_coordinates = raw_type(prev_active_set_);
// 		ctx->changed_coordinates     = raw_type(changed_);
		
// 		ctx->alpha = 1;
// 		ctx->blocksize = block_size_;
// 		// ctx->constrained = PETSC_TRUE;
// 		ctx->constrained = PETSC_FALSE;
// 		ctx->with_linesearch = PETSC_FALSE;


// 		if(block_size_ == 1) {
// 			ctx->constrain = &constrain_box;
// 		} else {
// 			ctx->constrain = &constrain_blocknd;
// 		}

// 		ctx_ = ctx;
// 	}
// }
