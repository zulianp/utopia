#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H

#include "utopia_petsc_Error.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Base.hpp"
#include "utopia_Size.hpp"
#include "utopia_Range.hpp"

#include "petscsys.h"
#include "petscmat.h"

#include <memory>

namespace utopia {
	class PetscMatrixMemory {
	public:
		PetscMatrixMemory(const MPI_Comm comm = PETSC_COMM_WORLD)
		: owner_(true)
		{
			MatCreate(comm, &_mat);
		}
		
		PetscMatrixMemory(Mat &mat, const bool owner = false)
		: _mat(mat), owner_(owner) { }
		
		~PetscMatrixMemory()
		{
			destroy();
		}

		inline void destroy()
		{
			if(owner_) {
				MatDestroy(&_mat);
			}
		}
		
		inline Mat &implementation() {
			return _mat;
		}
		
		inline const Mat &implementation() const {
            assert(_mat != nullptr);
			return _mat;
		}
		
		//MAT_DO_NOT_COPY_VALUES or MAT_COPY_VALUES, cause it to copy the numerical values in the matrix MAT_SHARE_NONZERO_PATTERN
		inline void duplicate(PetscMatrixMemory &other, MatDuplicateOption opt = MAT_COPY_VALUES) const {
			// if(other.owner_) {
				// MatDestroy(&other._mat);
			// }

			other.destroy();
			
			PetscErrorHandler::Check(MatDuplicate(_mat, opt, &other._mat));
		}
		
		inline void convert(PetscMatrixMemory &other, MatType newtype) {
			//MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
			PetscErrorHandler::Check(MatConvert(_mat, newtype, MAT_INITIAL_MATRIX, &other._mat));
		}
		
		inline void convert(MatType newtype) {
			//MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
			PetscErrorHandler::Check(MatConvert(_mat, newtype, MAT_REUSE_MATRIX, &_mat));
		}
		
		PetscMatrixMemory(const PetscMatrixMemory &other)
		: owner_(true)
		{
			PetscErrorHandler::Check(MatCopy(other._mat, _mat, SAME_NONZERO_PATTERN));
		}
		
		inline MPI_Comm communicator() const {
			MPI_Comm comm = PetscObjectComm((PetscObject) implementation());
			assert(comm != MPI_COMM_NULL);
			return comm;
		}
		
		void set_owner(const bool owner)
		{
			owner_ = owner;
		}
		
	private:
		Mat _mat;
		bool owner_;
	};
	
	class PetscMatrix {
	public:
		Mat &implementation() {
			return wrapper_->implementation();
		}

        inline MatType type() const
        {
            MatType type;
            MatGetType(implementation(), &type);
            return type;
        }

        inline std::string name() const
        {
            const char *name;
            PetscObjectGetName((PetscObject)implementation(), &name);
            return name;
        }

        inline void set_name(const std::string &name)
        {
            PetscObjectSetName((PetscObject)implementation(), name.c_str());
        }
		
		void wrap(Mat &mat)
		{
			wrapper_ = std::make_shared<PetscMatrixMemory>(mat, false);
		}
		
		const Mat &implementation() const {
			return wrapper_->implementation();
		}
		
		PetscMatrix(const MPI_Comm comm = PETSC_COMM_WORLD) {
			using std::make_shared;
			wrapper_ = make_shared<PetscMatrixMemory>(comm);
		}
		
		PetscMatrix(const PetscMatrix &other) {
			using std::make_shared;
			wrapper_ = make_shared<PetscMatrixMemory>();
			other.wrapper_->duplicate(*wrapper_);
		}
		
		PetscMatrix &operator=(const PetscMatrix &other) {
			if(wrapper_ == other.wrapper_) return *this;
			
			wrapper_ = std::make_shared<PetscMatrixMemory>();
			other.wrapper_->duplicate(*wrapper_);
			return *this;
		}
		
		PetscMatrix &operator=(PetscMatrix &&other) {
			if(wrapper_ == other.wrapper_) return *this;
			
			this->wrapper_ = other.wrapper_;
			other.wrapper_ = nullptr;
			return *this;
		}
		
		inline void describe() const {
			MatView(wrapper_->implementation(), PETSC_VIEWER_STDOUT_(communicator()));
		}
		
		inline MPI_Comm communicator() const {
			return wrapper_->communicator();
		}
		
		inline Size size() const
		{
			PetscInt rows, cols;
			check_error( MatGetSize(implementation(), &rows, &cols) );
			return {rows, cols};
		}
		
		inline Size local_size() const
		{
			PetscInt rows, cols;
			check_error( MatGetLocalSize(implementation(), &rows, &cols) );
			return {rows, cols};
		}
		
		inline PetscScalar get(const PetscInt row, const PetscInt col) const
		{
			PetscScalar value;
			check_error( MatGetValues(implementation(), 1, &row, 1, &col, &value) );
			return value;
		}
		
		inline void set(const PetscInt row, const PetscInt col, PetscScalar value) {
			check_error( MatSetValues(implementation(), 1, &row, 1, &col, &value, INSERT_VALUES) );
		}

		inline void add(const PetscInt row, const PetscInt col, PetscScalar value) {
			check_error( MatSetValues(implementation(), 1, &row, 1, &col, &value, ADD_VALUES) );
		}
		
		void add_matrix(const std::vector<PetscInt> &rows,
						const std::vector<PetscInt> &cols,
						const std::vector<PetscScalar> &values);
		
		void set_matrix(const std::vector<PetscInt> &rows,
						const std::vector<PetscInt> &cols,
						const std::vector<PetscScalar> &values);
		
		inline void write_lock() { }
		
		inline void write_unlock() {
			check_error( MatAssemblyBegin(implementation(), MAT_FINAL_ASSEMBLY) );
			check_error( MatAssemblyEnd(implementation(),   MAT_FINAL_ASSEMBLY) );
		}
		
		inline void read_lock() { }
		inline void read_unlock() { }

        void dense_init(
            MPI_Comm comm,
            MatType dense_type,
            PetscInt rows_local,
            PetscInt cols_local,
            PetscInt rows_global,
            PetscInt cols_global);


        bool write(const std::string &path) const;
        bool write_matlab(const std::string &path) const;

        bool read(MPI_Comm comm, const std::string &path);
        void copy_from(Mat mat);
        void copy_to(Mat mat);
        void copy_to(Mat *mat);

        inline Range row_range() const
        {
            PetscInt r_begin, r_end;
            MatGetOwnershipRange(implementation(), &r_begin, &r_end);
            assert(Range(r_begin, r_end).valid());
            return Range(r_begin, r_end);
        }

        inline Range col_range() const
        {
            PetscInt r_begin, r_end;
            MatGetOwnershipRangeColumn(implementation(), &r_begin, &r_end);
            assert(Range(r_begin, r_end).valid());
            return Range(r_begin, r_end);
        }

        void transpose();
        void transpose(PetscMatrix &result) const;

        void clear();
        bool is_sparse() const;

        void select(const std::vector<PetscInt> &row_index,
                    const std::vector<PetscInt> &col_index,
                    PetscMatrix &result) const;

        void local_select(
                const Range &local_row_range,
                const Range &local_col_range,
                const Range &global_col_range,
                PetscMatrix &result) const;

        void select(
                const Range &global_row_range,
                const Range &global_col_range,
                PetscMatrix &result) const;

        inline PetscScalar trace() const
        {
        	PetscScalar ret;
        	check_error( MatGetTrace(implementation(), &ret) );
        	return ret;
        }

        PetscScalar sum() const;
        PetscScalar max() const;
        PetscScalar min() const;

        void row_sum(PetscVector &col) const;
        void row_max(PetscVector &col) const;
        void row_min(PetscVector &col) const;

        inline PetscReal norm2() const
        {
        	PetscReal val;
        	MatNorm(implementation(), NORM_FROBENIUS, &val);
        	return val;
        }

        bool is_mpi() const;
        bool is_nan_or_inf() const;


        inline void scale(const PetscScalar factor)
        {
        	check_error( MatScale(implementation(), factor) );
        }

        //petsc says that it is correct only for square matrices
        void get_diag(PetscVector &result) const;
        void get_diag(PetscMatrix &result) const;

        inline void diag_shift(const PetscScalar factor)
        {
        	check_error( MatShift(implementation(), factor) );
        }

        void dense_init_diag(MatType dense_type, const PetscVector &diag);
        void matij_init_diag(const PetscVector &diag);

        void dense_init_values(
        	MPI_Comm comm,
        	MatType dense_type,
        	PetscInt local_rows,
        	PetscInt local_cols,
        	PetscInt global_rows,
        	PetscInt global_cols,
        	PetscScalar value
        );
       
        void dense_init_identity(
        	MPI_Comm comm,
        	MatType dense_type,
        	PetscInt local_rows,
        	PetscInt local_cols,
        	PetscInt global_rows,
        	PetscInt global_cols,
        	PetscScalar scale_factor);

        void matij_init_identity(
        	MPI_Comm comm,
        	PetscInt local_rows,
        	PetscInt local_cols,
        	PetscInt global_rows,
        	PetscInt global_cols,
        	PetscScalar scale_factor);

        void matij_init(
        	MPI_Comm comm,
        	PetscInt rows_local,
        	PetscInt cols_local,
        	PetscInt rows_global,
        	PetscInt cols_global,
        	PetscInt d_nnz,
        	PetscInt o_nnz
        );
        
		inline void destroy()
		{
			wrapper_->destroy();
		}

		void inverse(PetscMatrix &result) const;

		void mult(const PetscVector &vec, PetscVector &result) const;
		void mult_t(const PetscVector &vec, PetscVector &result) const;
		void mult(const PetscMatrix &mat, PetscMatrix &result) const;
		void mult_t(const PetscMatrix &mat, PetscMatrix &result) const;
		void mult_mat_t(const PetscMatrix &mat, PetscMatrix &result) const;

		///result = v2 + this * v1
		void mult_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const;
		void mult_t_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const;


		///this is y
		void axpy(const PetscScalar alpha, const PetscMatrix &x);

	private:
		std::shared_ptr<PetscMatrixMemory> wrapper_;
		
		inline static bool check_error(const PetscInt err)
		{
			return PetscErrorHandler::Check(err);
		}

        void select_aux(const std::vector<PetscInt> &row_index,
                        const std::vector<PetscInt> &col_index,
                        PetscMatrix &result) const;

        void par_assign_from_local_is(
                const std::vector<PetscInt> &remote_rows,
                const std::vector<PetscInt> &remote_cols,
                const PetscInt global_col_offset,
                const Range &local_col_range,
                PetscMatrix &result) const;

      
	};
}

#endif //UTOPIA_UTOPIA_PETSCMATRIX_H
