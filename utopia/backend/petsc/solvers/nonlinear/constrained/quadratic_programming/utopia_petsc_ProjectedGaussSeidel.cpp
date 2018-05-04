// #include "utopia_petsc_ProjectedGaussSeidel.hpp"

// #ifdef WITH_PASSO_EXTENSIONS

// #include "utopia_petsc_external.hpp"

// extern PetscErrorCode constrain_box(BLOCKSOLVER * blocksolver, MatScalar * diag, PetscScalar * t, PetscScalar * x,
// 	PetscScalar * lb, PetscScalar * ub, int row,
// 	int bs, PetscBool * constrained);

// extern PetscErrorCode constrain_blocknd(BLOCKSOLVER * blocksolver, MatScalar * diag, PetscScalar * t, PetscScalar * x,
// 	PetscScalar * lb, PetscScalar * ub, int row,
// 	int bs, PetscBool * constrained);

// namespace utopia {
	
// 	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::apply(const Vector &b, Vector &x)
// 	{
// 		if(this->verbose())
// 			this->init_solver("utopia/passo ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});

// 		const auto &A = *get_operator();

// 		auto ctx = std::static_pointer_cast<NBGS_CTX>(ctx_);

// 		if(empty(x) || size(x).get(0) != size(b).get(0)) {
// 			x = local_zeros(local_size(b));
// 		}

// 		old_x_ = x;
// 		bool converged = false;
// 		for(SizeType i = 0; i < this->max_it(); ++i) {
// 			auto err = NBGSStep(
// 					ctx.get(),
// 					raw_type(A),
// 					raw_type(b),
// 					raw_type(x),
// 					raw_type(*constraints_.lower_bound()),
// 					raw_type(*constraints_.upper_bound()),
// 					ctx->blocksize
// 				);

// 			if(err) return false;

// 			// std::cout << "na: " << ctx->num_active_coordinates << std::endl;

// 			PetscScalar diff = norm2(old_x_ - x);
// 			if(this->verbose())
// 				PrintInfo::print_iter_status({static_cast<Scalar>(i), diff}); 

// 			if(this->check_convergence(i, 1, 1, diff)) {
// 				converged = true;
// 				break;
// 			}

// 			prev_active_set_ = active_set_;
// 			old_x_ = x;
// 		}

// 		return converged;
// 	}

// 	bool ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::smooth(const Matrix &A, const Vector &b, Vector &x)
// 	{
// 		if(empty(x) || size(x).get(0) != size(b).get(0)) {
// 			x = local_zeros(local_size(b));
// 		}

// 		init(A);

// 		const Vector &ub = *constraints_.upper_bound();
// 		const Vector &lb = *constraints_.lower_bound();
// 		auto ctx = std::static_pointer_cast<NBGS_CTX>(ctx_);

// 		for(SizeType i = 0; i < this->sweeps(); ++i) {
	
// 			auto err = NBGSStep(
// 					ctx.get(),
// 					raw_type(A),
// 					raw_type(b),
// 					raw_type(x),
// 					raw_type(*constraints_.lower_bound()),
// 					raw_type(*constraints_.upper_bound()),
// 					ctx->blocksize
// 				);
// 		}

// 		return true;
// 	}

// 	ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::~ProjectedGaussSeidel() {}

// 	void ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC>::init(const Matrix &A)
// 	{
// 		auto ls = local_size(A);
// 		auto gs = size(A);

// 		DVectord dummy = local_zeros(ls.get(0));

// 		const bool has_bound = constraints_.has_bound();

// 		if(!constraints_.has_lower_bound()) {
// 			constraints_.lower_bound() = std::make_shared<Vector>(local_values(ls.get(0), -std::numeric_limits<PetscScalar>::max()));
// 		}

// 		if(!constraints_.has_upper_bound()) {
// 			constraints_.upper_bound() = std::make_shared<Vector>(local_values(ls.get(0), std::numeric_limits<PetscScalar>::max()));
// 		}

// 		auto ctx = std::shared_ptr<NBGS_CTX>(new NBGS_CTX, [](NBGS_CTX *&nbgs) { NBGSDestroy(nbgs); delete nbgs; nbgs = nullptr; });

// 		auto err = NBGSCreate(
// 			ctx.get(),
// 			raw_type(A),
// 			raw_type(dummy),
// 			raw_type(*constraints_.lower_bound()),
// 			raw_type(*constraints_.upper_bound()),
// 			block_size_,
// 			use_line_search_? PETSC_TRUE : PETSC_FALSE,
// 			has_bound? 		  PETSC_TRUE : PETSC_FALSE, 
// 			(block_size_ == 1)? &constrain_box :  &constrain_blocknd
// 		);

// 		ctx_ = ctx;
// 	}
// }

// #endif //WITH_PASSO_EXTENSIONS
