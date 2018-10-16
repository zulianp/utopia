#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"

#include <algorithm>
#include <set>

//PetscObjectTypeCompare((PetscObject)mat,newtype,&sametype);

namespace utopia {

    MatType PetscMatrix::type_override() const
    {
        return MATDENSE;
    }

    void PetscMatrix::add_matrix(
       const std::vector<PetscInt> &rows,
       const std::vector<PetscInt> &cols,
       const std::vector<PetscScalar> &values)
    {
        assert(rows.size() * cols.size() == values.size());

        check_error(
            MatSetValues(
               implementation(),
               static_cast<PetscInt>(rows.size()), &rows[0],
               static_cast<PetscInt>(cols.size()), &cols[0],
               &values[0],
               ADD_VALUES)
            );
    }

    void PetscMatrix::set_matrix(
       const std::vector<PetscInt> &rows,
       const std::vector<PetscInt> &cols,
       const std::vector<PetscScalar> &values)
    {
        assert(rows.size() * cols.size() == values.size());

        check_error(
            MatSetValues(
               implementation(),
               static_cast<PetscInt>(rows.size()), &rows[0],
               static_cast<PetscInt>(cols.size()), &cols[0],
               &values[0],
               INSERT_VALUES)
            );
    }

    void PetscMatrix::dense_init(
       MPI_Comm comm,
       MatType dense_type,
       PetscInt rows_local,
       PetscInt cols_local,
       PetscInt rows_global,
       PetscInt cols_global)
    {

        const std::string type_copy = dense_type;
        destroy();

        check_error( MatCreate(comm, &implementation()) );
        check_error( MatSetFromOptions(implementation()) );
        check_error( MatSetType(implementation(), type_copy.c_str()) );
        check_error( MatSetSizes(implementation(), rows_local, cols_local, rows_global, cols_global) );
        check_error( MatSetUp(implementation()) );

        check_error( MatSetOption(implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE) );
        check_error( MatSetOption(implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE) );
    }



    bool PetscMatrix::read(MPI_Comm comm, const std::string &path)
    {
        destroy();

        PetscViewer fd;

        bool err = check_error( PetscViewerBinaryOpen(comm, path.c_str(), FILE_MODE_READ, &fd) );
        err = err && check_error( MatCreate(comm, &implementation()) );
        err = err && check_error( MatSetType(implementation(), type_override()) );
        err = err && check_error( MatLoad(implementation(), fd) );

        check_error( PetscViewerDestroy(&fd) );
        return err;
    }

    bool PetscMatrix::write(const std::string &path) const
    {
        PetscViewer fd;

        bool err = check_error( PetscViewerBinaryOpen(communicator(), path.c_str(), FILE_MODE_WRITE, &fd) );
        err = err && check_error( MatView(implementation(), fd));

        check_error( PetscViewerDestroy(&fd) );
        return err;
    }

    bool PetscMatrix::write_matlab(const std::string &path) const
    {
        PetscViewer fd;

        bool err = check_error( PetscViewerASCIIOpen(communicator(), path.c_str(), &fd) );
        err = err && check_error( PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB) );
        err = err && check_error( MatView(implementation(), fd) );

        check_error( PetscViewerDestroy(&fd) );
        return err;
    }

    void PetscMatrix::copy_from(Mat mat)
    {
        destroy();
        check_error( MatDuplicate(mat, MAT_COPY_VALUES, &implementation()) );
    }

    void PetscMatrix::copy_to(Mat mat)
    {
        check_error( MatCopy(implementation(), mat, DIFFERENT_NONZERO_PATTERN) );
    }

    void PetscMatrix::copy_to(Mat *mat)
    {
        check_error( MatDuplicate(implementation(), MAT_COPY_VALUES, mat) );
    }

    void PetscMatrix::transpose()
    {
        check_error( MatTranspose(implementation(),  MAT_INPLACE_MATRIX, &implementation()) );
    }

    void PetscMatrix::transpose(PetscMatrix &result) const
    {
        if(implementation() == result.implementation()) {
            auto s = size();
            
            if(s.get(0) == s.get(1)) {
                result.transpose();
            } else {
                PetscMatrix temp;
                temp.destroy();

                check_error( MatTranspose(implementation(), MAT_INITIAL_MATRIX, &temp.implementation()) );
                result = std::move(temp);
            }

            return;
        }

        result.destroy();
        check_error( MatTranspose(implementation(), MAT_INITIAL_MATRIX, &result.implementation()) );
    }

    void PetscMatrix::clear()
    {
        MPI_Comm comm = communicator();
        destroy();
        MatCreate(comm, &implementation());
    }

    bool PetscMatrix::is_sparse() const
    {
        const std::string type_str(type());
        const size_t start = type_str.size() - 5;
        return !(type_str.substr(start, 5) == "dense");
    }

    void PetscMatrix::select(
       const std::vector<PetscInt> &row_index,
       const std::vector<PetscInt> &col_index,
       PetscMatrix &result) const
    {
        if(col_index.empty()) {
            PetscInt n_rows, n_cols;
            MatGetSize(implementation(), &n_rows, &n_cols);
            std::vector<PetscInt> all_index(n_cols);
            for(PetscInt i = 0; i < n_cols; ++i) {
                all_index[i] = i;
            }

            select_aux(row_index, all_index, result);

        } else {
            select_aux(row_index, col_index, result);
        }
    }

    void PetscMatrix::select_aux(const std::vector<PetscInt> &row_index,
       const std::vector<PetscInt> &col_index,
       PetscMatrix &result) const
    {
        // Mat r = implementation();
        MPI_Comm comm = communicator();

        PetscInt min_col = col_index[0], max_col = col_index[1];
        for(std::size_t i = 0; i < col_index.size(); ++i) {
            min_col = std::min(min_col, col_index[i]);
            max_col = std::max(max_col, col_index[i]);
        }

        int vals[2] = { min_col, -max_col };
        MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_INT, MPI_MIN, comm);
        min_col =  vals[0];
        max_col = -vals[1];

        PetscInt global_cols = max_col - min_col + 1;
        PetscInt local_cols = PETSC_DECIDE;
        PetscSplitOwnership(comm, &local_cols, &global_cols);

        unsigned long offsets_in = row_index.size();
        unsigned long offset_out = 0;

        MPI_Exscan(
         &offsets_in,
         &offset_out,
         1,
         MPI_UNSIGNED_LONG ,
         MPI_SUM,
         comm);

        offset_out += min_col;

        par_assign_from_local_is(
           row_index,
           col_index,
           min_col,
           Range(offset_out, offset_out + local_cols),
           result);
    }

    void PetscMatrix::par_assign_from_local_is(const std::vector<PetscInt> &remote_rows,
     const std::vector<PetscInt> &remote_cols,
     const PetscInt global_col_offset,
     const Range &local_col_range,
     PetscMatrix &result) const
    {
        MPI_Comm comm = communicator();
        Mat &l = result.implementation();
        const Mat r = implementation();

        int size;
        MPI_Comm_size(comm, &size);

        int rank;
        MPI_Comm_rank(comm, &rank);


        std::stringstream ss;

        for(auto r : remote_rows) {
            ss << r;
        }

        IS isrow;
        PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
        ierr = ISCreateGeneral(comm, remote_rows.size(), &remote_rows[0], PETSC_USE_POINTER, &isrow);

        IS iscol;
        ierr = ISCreateGeneral(comm, remote_cols.size(), &remote_cols[0], PETSC_USE_POINTER, &iscol);

        //TODO maybe make it with collective comms for finiding out if there are off proc entries
        bool has_off_proc_entries = size > 1;

        if(has_off_proc_entries) {
            Mat * l_ptr;

#if UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)
            ierr = MatGetSubMatrices(r, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, &l_ptr);
#else
            ierr = MatCreateSubMatrices(r, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, &l_ptr);
#endif

            MatType type;
            MatGetType(r, &type);
            MatSetType(l, type);
            MatSetSizes(l, remote_rows.size(), local_col_range.extent(), PETSC_DETERMINE, PETSC_DETERMINE);
            MatSetUp(l);

            PetscInt rbegin, rend, cbegin, cend;
            MatGetOwnershipRange(l, &rbegin, &rend);
            MatGetOwnershipRangeColumn(l, &cbegin, &cend);
            PetscInt n_rows = rend - rbegin;

            PetscInt n_values;
            const PetscInt * cols;
            const PetscScalar * values;


            for(PetscInt row = 0; row < n_rows; ++row) {

                MatGetRow(*l_ptr, row, &n_values, &cols, &values);

                for(PetscInt i = 0; i < n_values; ++i) {
                    MatSetValue(l, rbegin + row, cols[i] - global_col_offset, values[i], INSERT_VALUES);
                }

                MatRestoreRow(*l_ptr, row, &n_values, &cols, &values);
            }

            MatAssemblyBegin(l, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(l, MAT_FINAL_ASSEMBLY);

            MatDestroy(l_ptr);

        } else {
            result.destroy();
#if UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)
            ierr = MatGetSubMatrix(r, isrow, iscol, MAT_INITIAL_MATRIX, &l);
#else
            ierr = MatCreateSubMatrix(r, isrow, iscol, MAT_INITIAL_MATRIX, &l);
#endif //UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)
        }

        ISDestroy(&isrow);
        ISDestroy(&iscol);
    }

    void PetscMatrix::local_select(const Range &local_row_range,
     const Range &local_col_range,
     const Range &global_col_range,
     PetscMatrix &result) const
    {
        // PetscErrorCode ierr = 0;

        // Mat &l = result.implementation();
        // const Mat r = implementation();

        std::vector<PetscInt> remote_rows;
        remote_rows.reserve(local_row_range.extent());
        for(PetscInt l_row = local_row_range.begin(); l_row < local_row_range.end(); ++l_row) {
            remote_rows.push_back(l_row);
        }

        std::vector<PetscInt> remote_cols;
        remote_cols.reserve(global_col_range.extent());
        for(PetscInt l_col = global_col_range.begin(); l_col < global_col_range.end(); ++l_col) {
            remote_cols.push_back(l_col);
        }

        par_assign_from_local_is(remote_rows, remote_cols, global_col_range.begin(), local_col_range, result);
    }

    void PetscMatrix::select(const Range &global_row_range,
       const Range &global_col_range,
       PetscMatrix &result) const
    {

        const Mat r = implementation();
        PetscInt global_rows = global_row_range.extent();
        PetscInt local_rows  = PETSC_DECIDE;

        PetscInt global_cols = global_col_range.extent();
        PetscInt local_cols  = PETSC_DECIDE;

        MPI_Comm comm = PetscObjectComm((PetscObject)r);
        PetscSplitOwnership(comm, &local_rows, &global_rows);
        PetscSplitOwnership(comm, &local_cols, &global_cols);

        unsigned long offsets_in[2] = {
            (unsigned long)local_rows,
            (unsigned long)local_cols
        };

        unsigned long offset_out[2] = {
            (unsigned long)0,
            (unsigned long)0
        };

        MPI_Exscan(
         &offsets_in,
         &offset_out,
         2,
         MPI_UNSIGNED_LONG ,
         MPI_SUM,
         comm);

        offset_out[0] += global_row_range.begin();
        offset_out[1] += global_col_range.begin();

        local_select(
           Range(offset_out[0], offset_out[0] + local_rows),
           Range(offset_out[1], offset_out[1] + local_cols),
           global_col_range,
           result);
    }

    PetscScalar PetscMatrix::sum() const
    {
        Vec row_sum;

        MatCreateVecs(implementation(), nullptr, &row_sum);

        MatGetRowSum(implementation(), row_sum);

        PetscScalar res = 0.;
        check_error( VecSum(row_sum, &res) );

        VecDestroy(&row_sum);
        return res;
    }



    template<class Operation>
    inline static PetscScalar generic_local_reduce(const PetscMatrix &m, const PetscScalar &init_value, const Operation &op)
    {
        PetscScalar x = init_value;
        const PetscScalar * values;
        const PetscInt * cols;

        PetscInt r_begin, r_end;
        PetscInt n_values = 0;

        PetscInt local_r, local_c;

        MatGetLocalSize(m.implementation(), &local_r, &local_c);
        MatGetOwnershipRange(m.implementation(), &r_begin, &r_end);

        for(PetscInt row = r_begin; row < r_end; ++row) {

            MatGetRow(m.implementation(), row, &n_values, &cols, &values);

            if(n_values < local_c) {
                x = op.template apply<PetscScalar>(x, 0.);
            }

            for(PetscInt i = 0; i < n_values; ++i) {
                x = op.template apply<PetscScalar>(x, values[i]);
            }

            MatRestoreRow(m.implementation(), row, &n_values, &cols, &values);
        }

        return x;
    }

    PetscScalar PetscMatrix::max() const
    {
        PetscScalar result = -std::numeric_limits<PetscScalar>::max();
        MPI_Comm comm = communicator();

        int size = 0;
        MPI_Comm_size(comm, &size);

        if(size == 1 || !is_mpi()) {
            Vec v;
            MatCreateVecs(implementation(), nullptr, &v);
            MatGetRowMax(implementation(), v, nullptr);
            VecMax(v, nullptr, &result);

            VecDestroy(&v);
        } else {
            result = generic_local_reduce(*this, result, utopia::Max());
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX, comm);
        }

        return result;
    }

    PetscScalar PetscMatrix::min() const
    {
        PetscScalar result = std::numeric_limits<PetscScalar>::max();
        MPI_Comm comm = communicator();

        int size = 0;
        MPI_Comm_size(comm, &size);

        if(size == 1 || !is_mpi()) {
            Vec v;
            MatCreateVecs(implementation(), nullptr, &v);
            MatGetRowMin(implementation(), v, nullptr);
            VecMin(v, nullptr, &result);

            VecDestroy(&v);
        } else {
            result = generic_local_reduce(*this, result, utopia::Min());
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MIN, comm);
        }

        return result;
    }

    bool PetscMatrix::is_mpi() const
    {
        static const std::string seq = "seq";
        const std::string str = type();
        return std::search(begin(str), end(str), begin(seq), end(seq)) == str.end();
    }

    void PetscMatrix::get_diag(PetscVector &result) const
    {
        auto gs = size();

        result.destroy();

        if(gs.get(0) < gs.get(1)) {
            MatCreateVecs(implementation(), nullptr, &result.implementation());
        } else {
            MatCreateVecs(implementation(), &result.implementation(), nullptr);
        }

        check_error( MatGetDiagonal(implementation(), result.implementation()) );
        result.set_initialized(true);
    }

    void PetscMatrix::get_col(PetscVector &result, const PetscInt id) const
    {        
        auto gs = size();

        result.destroy();
        MatCreateVecs(implementation(), nullptr, &result.implementation());

        check_error( MatGetColumnVector(implementation(), result.implementation(), id) );
        result.set_initialized(true);
    }


    void PetscMatrix::dense_init_diag(MatType dense_type, const PetscVector &diag)
    {
        MPI_Comm comm = diag.communicator();
        PetscInt local_size  = diag.local_size();
        PetscInt global_size = diag.size();

        destroy();

        dense_init(
         comm,
         dense_type,
         local_size,
         local_size,
         global_size,
         global_size
         );

        check_error( MatZeroEntries(implementation()) );
        check_error( MatDiagonalSet( implementation(), diag.implementation(), INSERT_VALUES) );
    }

    void PetscMatrix::matij_init_diag(const PetscVector &diag)
    {
        auto local_size  = diag.local_size();
        auto global_size = diag.size();

        destroy();

        matij_init(
         diag.communicator(),
         type_override(),
         local_size,
         local_size,
         global_size,
         global_size,
         1,
         0);

        check_error( MatZeroEntries(implementation()) );
        check_error( MatDiagonalSet( implementation(), diag.implementation(), INSERT_VALUES) );
    }

    void PetscMatrix::nest(
       MPI_Comm comm,
       PetscInt nr,
       const IS is_row[],
       PetscInt nc,
       const IS is_col[],
       const Mat a[],
       const bool use_mat_nest_type
    )
    {
        destroy();

        if(use_mat_nest_type) {
            check_error( MatCreateNest(comm, nr, is_row, nc, is_col, a, &implementation()) );
        } else {
            Mat temp = nullptr;
            
            check_error( MatCreateNest(comm, nr, is_row, nc, is_col, a, &temp) );
            check_error( MatConvert(temp, type_override(), MAT_INITIAL_MATRIX, &implementation()) );

            check_error(  MatDestroy(&temp) );
        }
    }

    void PetscMatrix::get_diag(PetscMatrix &result) const
    {
        MatType type = this->type();
        MPI_Comm comm = communicator();

        result.destroy();
        check_error( MatCreate(comm, &result.implementation()) );
        check_error( MatSetType(result.implementation(), type) );

        const Size gs = size();
        const Size ls = local_size();
        // const bool is_row = gs.get(0) < gs.get(1);

        PetscVector vec;
        get_diag(vec);

        const PetscInt local_size = vec.local_size();
        const PetscInt global_size = vec.size();

        check_error( MatSetSizes(result.implementation(), local_size, local_size, global_size, global_size) );

        //in case it is a sparse format.

        //FIXME handle other cases
        check_error( MatSeqAIJSetPreallocation(result.implementation(), 1, PETSC_NULL) );
        check_error( MatMPIAIJSetPreallocation(result.implementation(), 1, PETSC_NULL, 0, PETSC_NULL) );
        check_error( MatSetUp(result.implementation()) );
        check_error( MatZeroEntries(result.implementation()) );

        check_error( MatDiagonalSet( result.implementation(), vec.implementation(), INSERT_VALUES ) );
    }

    bool PetscMatrix::empty() const
    {
        Size s = size();
        if(s.get(0) <= 0) return true;
        return false;
    }

    bool PetscMatrix::is_initialized_as(MPI_Comm comm,
       MatType dense_type,
       PetscInt local_rows,
       PetscInt local_cols,
       PetscInt global_rows,
       PetscInt global_cols)
    {
        if(empty()) {
            return false;
        }

        // TODO:: check type and comm

        PetscBool initialized;
        MatAssembled(implementation(), &initialized);

        if(initialized && (local_rows > 0 && global_cols > 0))
        {
            PetscInt m, n;
            MatGetLocalSize(implementation(), &m, &n);
            initialized = (m==local_rows && n == local_cols) ? PETSC_TRUE : PETSC_FALSE;
        }

        if(initialized)
        {
            PetscInt m, n;
            MatGetSize(implementation(), &m, &n);
            initialized = (m==global_rows && n == global_cols) ? PETSC_TRUE : PETSC_FALSE;
        }

        return initialized == PETSC_TRUE;
    }

    void PetscMatrix::dense_init_values(MPI_Comm comm,
        MatType dense_type,
        PetscInt local_rows,
        PetscInt local_cols,
        PetscInt global_rows,
        PetscInt global_cols,
        PetscScalar value
        )
    {
        if(!is_initialized_as(comm, dense_type, local_rows, local_cols, global_rows, global_cols)) {
            dense_init(comm, dense_type, local_rows, local_cols, global_rows, global_cols);
        }

        const auto r = row_range();
        const PetscInt r_begin = r.begin();
        const PetscInt r_end   = r.end();

        const PetscInt computed_global_cols = (global_cols <= 0) ? size().get(1) : global_cols;

        for (PetscInt i = r_begin; i < r_end; ++i) {
            for (PetscInt j = 0; j < computed_global_cols; ++j) {
                MatSetValue(implementation(), i, j, value, INSERT_VALUES);
            }
        }

        MatAssemblyBegin(implementation(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(implementation(), MAT_FINAL_ASSEMBLY);
    }


    void PetscMatrix::dense_init_identity(MPI_Comm comm,
      MatType dense_type,
      PetscInt local_rows,
      PetscInt local_cols,
      PetscInt global_rows,
      PetscInt global_cols,
      PetscScalar scale_factor)
    {
        if(!is_initialized_as(comm, dense_type, local_rows, local_cols, global_rows, global_cols))
            dense_init(comm, dense_type, local_rows, local_cols, global_rows, global_cols);

        check_error( MatZeroEntries(implementation()) );

        write_lock();

        const auto r = row_range();
        const PetscInt r_begin = r.begin();

        // otherwise global_cols gives -1, as it should be determined... 
        MatGetSize(implementation(), &global_rows, &global_cols); 
        const PetscInt r_end = PetscMin(r.end(), global_cols); 

        for(PetscInt i = r_begin; i < r_end; ++i) {
            set(i, i, scale_factor);
        }

        write_unlock();
    }

    // void PetscMatrix::matij_init_identity(
    //   MPI_Comm comm,
    //   PetscInt local_rows,
    //   PetscInt local_cols,
    //   PetscInt global_rows,
    //   PetscInt global_cols,
    //   PetscScalar scale_factor)
    // {
    //     matij_init_identity(comm, MATAIJ, local_rows, local_cols, global_rows, global_cols, scale_factor);
    // }

    void PetscMatrix::matij_init_identity(
      MPI_Comm comm,
      MatType sparse_type,
      PetscInt local_rows,
      PetscInt local_cols,
      PetscInt global_rows,
      PetscInt global_cols,
      PetscScalar scale_factor)
    {
        if(!is_initialized_as(comm, sparse_type, local_rows, local_cols, global_rows, global_cols))
            matij_init(
             comm,
             sparse_type,
             local_rows,
             local_cols,
             global_rows,
             global_cols,
             1,
             0
             );

        MatZeroEntries(implementation());

        write_lock();

        const auto r = row_range();
        const auto gs = size();

        const PetscInt r_begin = r.begin();
        const PetscInt r_end = PetscMin(r.end(), gs.get(1));

        for(PetscInt i = r_begin; i < r_end; ++i) {
            set(i, i, scale_factor);
        }

        write_unlock();

        // MatShift(implementation(), scale_factor);
    }

    // void PetscMatrix::matij_init(MPI_Comm comm,
    //    PetscInt rows_local,
    //    PetscInt cols_local,
    //    PetscInt rows_global,
    //    PetscInt cols_global,
    //    PetscInt d_nnz,
    //    PetscInt o_nnz)
    // {
    //     matij_init(comm, MATAIJ, rows_local, cols_local, rows_global, cols_global, d_nnz, o_nnz);
    // }

    void PetscMatrix::matij_init(MPI_Comm comm,
       MatType type,
       PetscInt rows_local,
       PetscInt cols_local,
       PetscInt rows_global,
       PetscInt cols_global,
       PetscInt d_nnz,
       PetscInt o_nnz)
    {
        destroy();

        check_error( MatCreate(comm, &implementation()) );
        check_error( MatSetSizes(implementation(), rows_local, cols_local, rows_global, cols_global) );

        check_error( MatSetType(implementation(), type) );
        check_error( MatSeqAIJSetPreallocation(implementation(), PetscMax(d_nnz, 1), PETSC_NULL) );
        check_error( MatMPIAIJSetPreallocation(implementation(), PetscMax(d_nnz, 1), PETSC_NULL, PetscMax(o_nnz, 1), PETSC_NULL) );

        check_error( MatSetOption(implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE) );
        check_error( MatSetOption(implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE) );
        check_error( MatSetOption(implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE) );

        check_error( MatZeroEntries(implementation()) );
    }


    // void PetscMatrix::matij_init(MPI_Comm comm,
    //    PetscInt rows_local,
    //    PetscInt cols_local,
    //    PetscInt rows_global,
    //    PetscInt cols_global,
    //    const std::vector<PetscInt> &d_nnz,
    //    const std::vector<PetscInt> &o_nnz)
    // {
    //     matij_init(comm, MATAIJ, rows_local, cols_local, rows_global, cols_global, d_nnz, o_nnz);
    // }


    void PetscMatrix::matij_init(MPI_Comm comm,
       MatType type,
       PetscInt rows_local,
       PetscInt cols_local,
       PetscInt rows_global,
       PetscInt cols_global,
       const std::vector<PetscInt> &d_nnz,
       const std::vector<PetscInt> &o_nnz)
    {
        destroy();

        check_error( MatCreate(comm, &implementation()) );
        check_error( MatSetSizes(implementation(), rows_local, cols_local, rows_global, cols_global) );

        check_error( MatSetType(implementation(), type) );
        check_error( MatSeqAIJSetPreallocation(implementation(), PETSC_DEFAULT , &d_nnz[0]) );
        check_error( MatMPIAIJSetPreallocation(implementation(), PETSC_DEFAULT , &d_nnz[0], PETSC_DEFAULT, &o_nnz[0]) );

        check_error( MatSetOption(implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE) );
        check_error( MatSetOption(implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE) );
        check_error( MatSetOption(implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE) );

        check_error( MatZeroEntries(implementation()) );
    }

    // void PetscMatrix::mat_aij_cusparse_init(MPI_Comm comm,
    //     PetscInt rows_local,
    //     PetscInt cols_local,
    //     PetscInt rows_global,
    //     PetscInt cols_global,
    //     PetscInt d_nnz,
    //     PetscInt o_nnz
    //     )
    // {
    //     matij_init(comm, MATAIJCUSPARSE, rows_local, cols_local, rows_global, cols_global, d_nnz, o_nnz);
    // }

    void PetscMatrix::mat_baij_init(MPI_Comm comm,
        PetscInt rows_local,
        PetscInt cols_local,
        PetscInt rows_global,
        PetscInt cols_global,
        PetscInt d_nnz,
        PetscInt o_nnz,
        PetscInt block_size
        )
    {
        destroy();

        check_error( MatCreate(comm, &implementation()) );
        check_error( MatSetSizes(implementation(), rows_local, cols_local, rows_global, cols_global) );
        check_error( MatSetBlockSize(implementation(), block_size));

        check_error( MatSetType(implementation(), MATBAIJ) );
        check_error( MatSeqBAIJSetPreallocation(implementation(), block_size, PetscMax(d_nnz, 1), PETSC_NULL) );
        check_error( MatMPIBAIJSetPreallocation(implementation(), block_size, PetscMax(d_nnz, 1), PETSC_NULL, PetscMax(o_nnz, 1), PETSC_NULL) );

        check_error( MatSetOption(implementation(), MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE) );
        check_error( MatSetOption(implementation(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE) );
        check_error( MatSetOption(implementation(), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE) );

        check_error( MatZeroEntries(implementation()) );
    }

    bool PetscMatrix::is_nan_or_inf() const
    {
        int has_nan = 0;
        const PetscScalar * values;
        const PetscInt * cols;

        PetscInt r_begin, r_end;
        PetscInt n_values = 0;

        PetscInt local_r, local_c;
        MatGetLocalSize(implementation(), &local_r, &local_c);
        MatGetOwnershipRange(implementation(), &r_begin, &r_end);

        for(PetscInt row = r_begin; row < r_end; ++row) {

            MatGetRow(implementation(), row, &n_values, &cols, &values);

            for(PetscInt i = 0; i < n_values; ++i) {
                has_nan = (PetscIsInfOrNanScalar(values[i])? 1 : 0);
                if(has_nan) break;
            }

            MatRestoreRow(implementation(), row, &n_values, &cols, &values);
            if(has_nan) break;
        }

        if(is_mpi()) {
            MPI_Allreduce(MPI_IN_PLACE, &has_nan, 1, MPI_INT, MPI_MAX, communicator());
        }

        return has_nan > 0;
    }

    void PetscMatrix::inverse(PetscMatrix &result) const
    {
        assert(!is_sparse());

        PetscMatrix I, L = *this;

        Size gs = size();
        Size ls = local_size();

        I.dense_init(
           communicator(),
           type(),
           ls.get(0),
           ls.get(1),
           gs.get(0),
           gs.get(1)
           );

        check_error( MatZeroEntries(I.implementation()) );
        check_error( MatShift(I.implementation(), 1.) );

        result.dense_init(
          communicator(),
          type(),
          ls.get(0),
          ls.get(1),
          gs.get(0),
          gs.get(1)
          );

        check_error( MatZeroEntries(result.implementation()) );

        IS isr, isc;
        MatFactorInfo info;

        check_error( MatGetOrdering(L.implementation(), MATORDERINGNATURAL, &isr, &isc) );
        check_error( MatLUFactor( L.implementation(), isr, isc, &info ) );
        check_error( MatMatSolve(L.implementation(), I.implementation(), result.implementation()) );


        check_error( ISDestroy(&isr) );
        check_error( ISDestroy(&isc) );
    }

    template<class Operation>
    inline static void reduce_rows(PetscVector &result,
     const PetscMatrix &mat,
     const PetscScalar &init_value,
     const Operation &op
     )
    {
        assert(!result.is_null());

        const PetscScalar * values;
        const PetscInt * cols;

        PetscInt r_begin, r_end;
        PetscInt n_values = 0;

        PetscInt global_r, global_c, local_r, local_c;
        MatGetSize(mat.implementation(), &global_r, &global_c);
        MatGetLocalSize(mat.implementation(), &local_r, &local_c);

        MatGetOwnershipRange(mat.implementation(), &r_begin, &r_end);

        result.write_lock();

        for(PetscInt row = r_begin; row < r_end; ++row) {
            MatGetRow(mat.implementation(), row, &n_values, &cols, &values);

            PetscScalar x = init_value;
            for(PetscInt i = 0; i < n_values; ++i) {
                x = op.template apply<PetscScalar>(x, values[i]);
            }

            if(n_values < global_c) {
                x = op.template apply<PetscScalar>(x, 0.);
            }

            MatRestoreRow(mat.implementation(), row, &n_values, &cols, &values);
            VecSetValues(result.implementation(), 1, &row, &x, INSERT_VALUES);
        }

        result.write_unlock();
    }

    void PetscMatrix::row_sum(PetscVector &col) const
    {
        MPI_Comm comm = communicator();

        if(col.is_null() || col.size() != size().get(0)) {
            col.destroy();
            MatCreateVecs(implementation(), nullptr, &col.implementation());
        }

        int size = 0;
        MPI_Comm_size(comm, &size);

        if(size == 1 || !is_mpi()) {
            MatGetRowSum(implementation(), col.implementation());
            col.set_initialized(true);
        } else {
            reduce_rows(
                col,
                *this,
                0.,
                utopia::Plus()
                );
        }
    }

    void PetscMatrix::row_max(PetscVector &col) const
    {
        MPI_Comm comm = communicator();

        if(col.is_null() || col.size() != size().get(0)) {
            col.destroy();
            MatCreateVecs(implementation(), nullptr, &col.implementation());
        }

        int size = 0;
        MPI_Comm_size(comm, &size);

        if(size == 1 || !is_mpi()) {
            MatGetRowMax(implementation(), col.implementation(), nullptr);
        } else {
            reduce_rows(
                col,
                *this,
                -std::numeric_limits<PetscScalar>::max(),
                utopia::Max()
                );
        }
    }

    void PetscMatrix::row_min(PetscVector &col) const
    {
        MPI_Comm comm = communicator();

        if(col.is_null() || col.size() != size().get(0)) {
            col.destroy();
            MatCreateVecs(implementation(), nullptr, &col.implementation());
        }

        int size = 0;
        MPI_Comm_size(comm, &size);

        if(size == 1 || !is_mpi()) {
            MatGetRowMin(implementation(), col.implementation(), nullptr);
        } else {
            reduce_rows(
                col,
                *this,
                std::numeric_limits<PetscScalar>::max(),
                utopia::Min()
                );
        }
    }

    void PetscMatrix::col_sum(PetscVector &col) const
    {
        PetscVector temp;
        temp.values(communicator(), col.type(), local_size().get(0), size().get(0), 1.);
        this->mult_t(temp, col);
    }

    void PetscMatrix::mult(const PetscVector &vec, PetscVector &result) const
    {
        if(vec.implementation() == result.implementation()) {
            assert(false && "handle me");
        }

        assert(vec.is_consistent());
        assert(result.is_consistent());

        MPI_Comm comm = vec.communicator();

        if(result.is_null()) {
            // MatCreateVecs(implementation(), nullptr, &result.implementation());
            create_vecs(nullptr, &result.implementation());
        } else if(comm != result.communicator()) {
            result.destroy();
            // MatCreateVecs(implementation(), nullptr, &result.implementation());
            create_vecs(nullptr, &result.implementation());
        } else {
            Size gs = size();
            Size ls = local_size();
            VecSetSizes(result.implementation(), ls.get(0), gs.get(0));
        }

        check_error( MatMult(implementation(), vec.implementation(), result.implementation() ) );

        result.set_initialized(true);

        assert(result.implementation() != nullptr);
        assert(result.is_consistent());
        // assert(result.same_type(vec)); //FIXME
    }

    void PetscMatrix::mult_t(const PetscVector &vec, PetscVector &result) const
    {
        if(vec.implementation() == result.implementation()) {
            assert(false && "handle me");
        }

        MPI_Comm comm = vec.communicator();

        if(result.is_null()) {
            MatCreateVecs(implementation(), &result.implementation(), nullptr);
        } else if(comm != result.communicator()) {
            result.destroy();
            MatCreateVecs(implementation(), &result.implementation(), nullptr);
        } else {
            Size gs = size();
            Size ls = local_size();
            VecSetSizes(result.implementation(), ls.get(1), gs.get(1));
        }

        check_error( MatMultTranspose(implementation(), vec.implementation(), result.implementation() ) );
        assert(result.implementation() != nullptr);
        result.set_initialized(true);
    }

    void PetscMatrix::mult(const PetscMatrix &mat, PetscMatrix &result) const
    {
        PetscBool      flg;
        // this is very unefficient hack, but still better than fail... 
        PetscObjectTypeCompareAny((PetscObject)mat.implementation(),&flg,MATMPIDENSE,NULL);
        if (flg)
        {
            if(mat.implementation() != result.implementation() && implementation() != result.implementation())
            {
                result.destroy();
                Mat temp; 
                MatConvert(mat.implementation(), MATMPIAIJ, MAT_INITIAL_MATRIX, &temp); 
                MatMatMult(implementation(), temp, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation());
                MatDestroy(&temp); 
            }
            else 
            {
                PetscMatrix temp2; 
                temp2.destroy(); 

                Mat temp; 
                MatConvert(mat.implementation(), MATMPIAIJ, MAT_INITIAL_MATRIX, &temp); 
                MatMatMult(implementation(), temp, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2.implementation());
                MatDestroy(&temp); 
                result = std::move(temp2);
            }                
        }
        else
        {
            if(mat.implementation() != result.implementation() && implementation() != result.implementation()) 
            {
                result.destroy();
                MatMatMult(implementation(), mat.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation());
            } else {
                PetscMatrix temp;
                temp.destroy();
                MatMatMult(implementation(), mat.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.implementation());
                result = std::move(temp);
            }
        }
    }

    void PetscMatrix::mult_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const
    {
        if (v1.implementation() == result.implementation() || v2.implementation() == result.implementation()) {
            PetscVector temp;
            temp.repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size());
            MatMultAdd(implementation(), v1.implementation(), v2.implementation(), temp.implementation());
            result = std::move(temp);
        } else {
            result.repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size());
            MatMultAdd(implementation(), v1.implementation(), v2.implementation(), result.implementation());
        }
    }

    void PetscMatrix::mult_t_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const
    {
        if (v1.implementation() == result.implementation() || v2.implementation() == result.implementation()) {
            PetscVector temp;
            temp.repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size());
            MatMultTransposeAdd(implementation(), v1.implementation(), v2.implementation(), temp.implementation());
            result = std::move(temp);
        } else {
            result.repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size());
            MatMultTransposeAdd(implementation(), v1.implementation(), v2.implementation(), result.implementation());
        }
    }

    void PetscMatrix::mult_t(const PetscMatrix &mat, PetscMatrix &result) const
    {
        if(mat.implementation() != result.implementation() && implementation() != result.implementation()) {
            result.destroy();
            check_error( MatTransposeMatMult(implementation(), mat.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) );
        } else {
            PetscMatrix temp; temp.destroy();
            check_error( MatTransposeMatMult(implementation(), mat.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.implementation()) );
            result = std::move(temp);
        }
    }

    void PetscMatrix::mult_mat_t(const PetscMatrix &mat, PetscMatrix &result) const
    {
        m_utopia_warning("> FIXME MatMatTransposeMult does not work in parallel, prepare work around");

        if(mat.implementation() != result.implementation() && implementation() != result.implementation()) {
            result.destroy();
            check_error( MatMatTransposeMult(implementation(), mat.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.implementation()) );
        } else {
            PetscMatrix temp; temp.destroy();
            check_error( MatMatTransposeMult(implementation(), mat.implementation(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.implementation()) );
            result = std::move(temp);
        }
    }

    void PetscMatrix::axpy(const PetscScalar alpha, const PetscMatrix &x) {
        check_error( MatAXPY(implementation(), alpha, x.implementation(), DIFFERENT_NONZERO_PATTERN) );
    }

    void PetscMatrix::convert_to_mat_baij(const PetscInt block_size)
    {
        auto ls = local_size();
        auto gs = size();

        PetscMatrix temp;
        temp.mat_baij_init(
         communicator(),
         ls.get(0),
         ls.get(1),
         gs.get(0),
         gs.get(1),
         PETSC_DEFAULT,
         PETSC_DEFAULT,
         block_size
         );

        check_error( MatCopy(implementation(), temp.implementation(), DIFFERENT_NONZERO_PATTERN) );

        *this = std::move(temp);
    }

    //testing MATAIJCUSPARSE,MATSEQAIJCUSPARSE
    bool PetscMatrix::PetscMatrix::is_cuda() const
    {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare((PetscObject) implementation(), MATAIJCUSPARSE, &match);
        if(match == PETSC_TRUE) return true;

        PetscObjectTypeCompare((PetscObject) implementation(), MATSEQAIJCUSPARSE, &match);
        return match == PETSC_TRUE;
    }

    VecType PetscMatrix::compatible_cuda_vec_type() const
    {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare((PetscObject) implementation(), MATAIJCUSPARSE, &match);

        if(match) {
            return VECMPICUDA;
        } else {
            return VECSEQCUDA;
        }
    }

    bool PetscMatrix::create_vecs(Vec *x, Vec *y) const
    {
        // if(is_cuda()) {
        //     auto vec_type = compatible_cuda_vec_type();
        //     Size s = size();
        //     Size ls = local_size();

        //     MPI_Comm comm = communicator();

        // 	if(x) {
        //         check_error( VecCreate(comm, x) );
        //         check_error( VecSetFromOptions(*x) );
        //         check_error( VecSetType(*x, vec_type) );

        //         check_error( VecSetSizes(*x, ls.get(1), s.get(1)) );
        // 	}

        // 	if(y) {
        //         check_error( VecCreate(comm, y) );
        //         check_error( VecSetFromOptions(*y) );
        //         check_error( VecSetType(*y, vec_type) );
        //         check_error( VecSetSizes(*y, ls.get(0), s.get(0)) );
        // 	}

        // } else {

        	MatCreateVecs(implementation(), x, y);
        // }

        return true;
    }

    bool PetscMatrix::has_type(VecType type) const
    {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare((PetscObject) implementation(), type, &match);
        return match == PETSC_TRUE;
    }

    bool PetscMatrix::same_type(const PetscMatrix &other) const
    {
        return has_type(other.type());
    }
}
