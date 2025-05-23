#include "utopia_petsc_Matrix.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Operators.hpp"
#include "utopia_TypeToString.hpp"
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_petsc_CrsView.hpp"

#include <algorithm>
#include <set>
#include <utility>

#ifdef UTOPIA_ENABLE_MATRIX_IO

#include "matrixio_crs.h"
#include "utils.h"

#endif

#ifdef UTOPIA_ENABLE_YAML_CPP
#include "utopia_YAMLInput.hpp"
#endif

// PetscObjectTypeCompare((PetscObject)mat,newtype,&sametype);
// Experts: Mark Hoemmen, Chris Siefert

namespace utopia {

    void PetscMatrix::copy(const PetscMatrix &other) /*override*/
    {
        if (raw_type() == other.raw_type()) {
            if (this == &other) {
                return;
            }

            // This could happen only when wrapping (i.e. using the wrap method)
            m_utopia_warning("should we just duplicate its memory?") return;
        }

        if (other.empty()) {
            clear();
            return;
        }

        wrapper_ = std::make_shared<PetscMatrixMemory>();
        other.wrapper_->duplicate(*wrapper_);
    }

    void PetscMatrix::same_nnz_pattern_copy(const PetscMatrix &other) {
        assert(!other.empty());
        MatCopy(other.raw_type(), raw_type(), SAME_NONZERO_PATTERN);
    }

    void PetscMatrix::transform(const Sqrt &op) { op_transform(op); }

    void PetscMatrix::transform(const Pow2 &op) { op_transform(op); }

    void PetscMatrix::transform(const Log &op) { op_transform(op); }

    void PetscMatrix::transform(const Exp &op) { op_transform(op); }

    void PetscMatrix::transform(const Cos &op) { op_transform(op); }

    void PetscMatrix::transform(const Sin &op) { op_transform(op); }

    void PetscMatrix::transform(const Abs &op) { op_transform(op); }

    void PetscMatrix::transform(const Minus &op) { op_transform(op); }

    void PetscMatrix::transform(const Pow &op) { op_transform(op); }

    void PetscMatrix::transform(const Reciprocal<Scalar> &op) { op_transform(op); }

    void PetscMatrix::transform(std::function<Scalar(const Scalar &)> op) { transform_values(std::move(op)); }

    void PetscMatrix::transform(std::function<Scalar(const SizeType &, const SizeType &, const Scalar &)> op) {
        transform_ijv(op);
    }

    MatType PetscMatrix::type_override() const {
        return MATAIJ;
        // return MATDENSE;
    }

    void PetscMatrix::add_matrix(const std::vector<PetscInt> &rows,
                                 const std::vector<PetscInt> &cols,
                                 const std::vector<Scalar> &values) {
        assert(rows.size() * cols.size() == values.size());

        check_error(MatSetValues(raw_type(),
                                 static_cast<PetscInt>(rows.size()),
                                 &rows[0],
                                 static_cast<PetscInt>(cols.size()),
                                 &cols[0],
                                 &values[0],
                                 ADD_VALUES));
    }

    void PetscMatrix::set_matrix(const std::vector<PetscInt> &rows,
                                 const std::vector<PetscInt> &cols,
                                 const std::vector<Scalar> &values) {
        assert(rows.size() * cols.size() == values.size());

        check_error(MatSetValues(raw_type(),
                                 static_cast<PetscInt>(rows.size()),
                                 &rows[0],
                                 static_cast<PetscInt>(cols.size()),
                                 &cols[0],
                                 &values[0],
                                 INSERT_VALUES));
    }

    void PetscMatrix::dense_init(MPI_Comm comm,
                                 MatType dense_type,
                                 PetscInt rows_local,
                                 PetscInt cols_local,
                                 PetscInt rows_global,
                                 PetscInt cols_global) {
        const std::string type_copy = dense_type;
        destroy();

        check_error(MatCreate(comm, &raw_type()));
        check_error(MatSetFromOptions(raw_type()));
        check_error(MatSetType(raw_type(), type_copy.c_str()));
        check_error(MatSetSizes(raw_type(), rows_local, cols_local, rows_global, cols_global));
        check_error(MatSetUp(raw_type()));

        check_error(MatSetOption(raw_type(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE));
        check_error(MatSetOption(raw_type(), MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE));

        UTOPIA_REPORT_ALLOC("PetscMatrix::dense_init");
    }

    bool PetscMatrix::read(MPI_Comm comm, const std::string &path) {
#ifdef UTOPIA_ENABLE_MATRIX_IO
        Path ppath = path;

#ifdef UTOPIA_ENABLE_YAML_CPP
        auto ext = ppath.extension();
        if (ext == "yaml" || ext == "yml") {
            Path folder = ppath.parent();
            Path rowptr_path = folder / "rowptr.raw";
            Path colidx_path = folder / "colidx.raw";
            Path values_path = folder / "values.raw";

            std::string rowptr_type = TypeToString<SizeType>::get();
            std::string colidx_type = TypeToString<SizeType>::get();
            std::string values_type = TypeToString<Scalar>::get();

            YAMLInput yaml;
            if (yaml.open(path)) {
                return false;
            }

            // Read paths from file
            yaml.get("rowptr_path", rowptr_path);
            yaml.get("colidx_path", colidx_path);
            yaml.get("values_path", values_path);

            // Read type information from file
            yaml.get("rowptr_type", rowptr_type);
            yaml.get("colidx_type", colidx_type);
            yaml.get("values_type", values_type);

            return read_raw(comm, rowptr_path, colidx_path, values_path, rowptr_type, colidx_type, values_type);
        }
#endif

        if (ppath.extension() == "raw") {
            Path folder = ppath.parent();
            return read_raw(comm, folder / "rowptr.raw", folder / "colidx.raw", folder / "values.raw");
        }
#endif  // UTOPIA_ENABLE_MATRIX_IO

        destroy();

        PetscViewer fd;

        bool err = check_error(PetscViewerBinaryOpen(comm, path.c_str(), FILE_MODE_READ, &fd));
        err = err && check_error(MatCreate(comm, &raw_type()));
        err = err && check_error(MatSetType(raw_type(), type_override()));
        err = err && check_error(MatLoad(raw_type(), fd));

        check_error(PetscViewerDestroy(&fd));

        UTOPIA_REPORT_ALLOC("PetscMatrix::read");
        return err;
    }

#ifdef UTOPIA_ENABLE_MATRIX_IO
    bool PetscMatrix::read_raw(MPI_Comm comm,
                               const std::string &rowptr_path,
                               const std::string &colidx_path,
                               const std::string &values_path,
                               const std::string &rowptr_type,
                               const std::string &colidx_type,
                               const std::string &values_type) {
        crs_t crs;
        if (crs_read_str(comm,
                         rowptr_path.c_str(),
                         colidx_path.c_str(),
                         values_path.c_str(),
                         rowptr_type.c_str(),
                         colidx_type.c_str(),
                         values_type.c_str(),
                         &crs)) {
            return false;
        }

        // if()

        PetscInt *rowptr = (PetscInt *)crs.rowptr;
        PetscInt firstrow = rowptr[0];

        // Remove global indexing
        for (ptrdiff_t i = 0; i <= crs.lrows; ++i) {
            rowptr[i] -= firstrow;
        }

        // FIXME add workaround
        // assert(crs.values_type_size == sizeof(Scalar));
        // assert(crs.rowptr_type_size == sizeof(SizeType));
        // assert(crs.colidx_type_size == sizeof(SizeType));

        destroy();

        int size;
        MPI_Comm_size(comm, &size);

        PetscInt lcols = 0;
        PetscInt *colidx = (PetscInt *)crs.colidx;

        for (ptrdiff_t i = 0; i < crs.lnnz; ++i) {
            lcols = std::max(colidx[i], lcols);
        }

        lcols += 1;

        if (size == 1) {
            check_error(MatCreateSeqAIJWithArrays(
                comm, crs.grows, lcols, rowptr, (PetscInt *)crs.colidx, (Scalar *)crs.values, &raw_type()));

        } else {
            PetscInt gcols = lcols;
            MPI_Allreduce(MPI_IN_PLACE, &gcols, 1, MPIType<PetscInt>::value(), MPI_MAX, comm);

            if (gcols == crs.grows) {
                lcols = (PetscInt)crs.lrows;
            } else {
                Utopia::Abort("read_raw does not support rectangular matrices!");
            }

            check_error(MatCreateMPIAIJWithArrays(comm,
                                                  crs.lrows,
                                                  lcols,
                                                  crs.grows,
                                                  gcols,
                                                  (PetscInt *)crs.rowptr,
                                                  (PetscInt *)crs.colidx,
                                                  (Scalar *)crs.values,
                                                  &raw_type()));
        }

        destroy_callback = [crs]() {
            crs_t crs_copy = crs;
            crs_free(&crs_copy);
        };

        UTOPIA_REPORT_ALLOC("PetscMatrix::read_raw");
        return true;
    }

    bool PetscMatrix::read_raw(MPI_Comm comm,
                               const std::string &rowptr_path,
                               const std::string &colidx_path,
                               const std::string &values_path) {
        return read_raw(comm,
                        rowptr_path,
                        colidx_path,
                        values_path,
                        TypeToString<SizeType>::get(),
                        TypeToString<SizeType>::get(),
                        TypeToString<Scalar>::get());
    }

    bool PetscMatrix::write_raw(MPI_Comm comm,
                                const std::string &rowptr_path,
                                const std::string &colidx_path,
                                const std::string &values_path) const {
        crs_t crs;
        crs.lrows = this->local_rows();
        crs.grows = this->rows();

        crs.rowptr_type = string_to_mpi_datatype(TypeToString<SizeType>::get());
        crs.colidx_type = string_to_mpi_datatype(TypeToString<SizeType>::get());
        crs.values_type = string_to_mpi_datatype(TypeToString<Scalar>::get());

        int size;
        MPI_Comm_size(comm, &size);
        if (size == 1) {
            PetscCrsView crs_view(raw_type());

            crs.rowptr = (char *)&crs_view.row_ptr()[0];
            crs.colidx = (char *)&crs_view.colidx()[0];
            crs.values = (char *)&crs_view.values()[0];

            crs.lnnz = crs_view.nnz();
            crs.gnnz = crs_view.nnz();

            crs.start = 0;
            crs.rowoffset = 0;

        } else {
            PetscErrorCode err = 0;

            const PetscInt *cols;
            Mat d, o;
            err = MatMPIAIJGetSeqAIJ(raw_type(), &d, &o, &cols);
            assert(err == 0);

            PetscCrsView d_crs_view(d);
            PetscCrsView o_crs_view(o);

            crs.lnnz = d_crs_view.nnz() + o_crs_view.nnz();
            PetscInt *rowptr = (PetscInt *)malloc((crs.lrows + 1) * sizeof(PetscInt));
            PetscInt *colidx = (PetscInt *)malloc(crs.lnnz * sizeof(PetscInt));
            PetscScalar *values = (PetscScalar *)malloc(crs.lnnz * sizeof(PetscScalar));

            long lnnz = crs.lnnz;
            long start = 0;
            MPI_Exscan(&lnnz, &start, 1, MPI_LONG, MPI_SUM, comm);
            crs.start = start;
            crs.rowoffset = this->row_range().begin();

            PetscInt coloff = this->col_range().begin();

            auto d_rowptr = d_crs_view.row_ptr();
            auto o_rowptr = o_crs_view.row_ptr();

            auto d_colidx = d_crs_view.colidx();
            auto o_colidx = o_crs_view.colidx();

            auto d_values = d_crs_view.values();
            auto o_values = o_crs_view.values();

            // PetscInt nghosts;
            // const PetscInt *ghosts = nullptr;
            // MatGetGhosts(raw_type(), &nghosts, &ghosts);

            rowptr[0] = start;

            for (ptrdiff_t i = 0; i < crs.lrows; i++) {
                rowptr[i + 1] = (d_rowptr[i + 1] - d_rowptr[i]) + (o_rowptr[i + 1] - o_rowptr[i]) + rowptr[i];

                ptrdiff_t last_offset = rowptr[i];
                for (PetscInt k = d_rowptr[i]; k < d_rowptr[i + 1]; k++) {
                    colidx[last_offset + k] = coloff + d_colidx[k];
                    values[last_offset + k] = d_values[k];

                    assert(colidx[last_offset] < this->cols());
                }

                last_offset += d_rowptr[i + 1] - d_rowptr[i];
                for (PetscInt k = o_rowptr[i]; k < o_rowptr[i + 1]; k++) {
                    colidx[last_offset + k] = cols[o_colidx[k]];
                    values[last_offset + k] = o_values[k];

                    assert(colidx[last_offset] < this->cols());
                }
            }

            crs.rowptr = (char *)rowptr;
            crs.colidx = (char *)colidx;
            crs.values = (char *)values;
        }

        bool ok = crs_write(comm, rowptr_path.c_str(), colidx_path.c_str(), values_path.c_str(), &crs) == 0;

        if (size == 1) {
            crs_release(&crs);
        } else {
            crs_free(&crs);
        }

        return ok;
    }

#endif  // UTOPIA_ENABLE_MATRIX_IO

    bool PetscMatrix::write(const std::string &path) const {
#ifdef UTOPIA_ENABLE_MATRIX_IO
        Path ppath = path;
        if (ppath.extension() == "raw") {
            Path folder = ppath.parent();
            folder.make_dir();
            return write_raw(comm().get(), folder / "rowptr.raw", folder / "colidx.raw", folder / "values.raw");
        }
#endif

        if (is_matlab_file(path)) {
            return write_matlab(path);
        }
        return write_binary(path);
    }

    bool PetscMatrix::write_binary(const std::string &path) const {
        PetscViewer fd;

        bool err = check_error(PetscViewerBinaryOpen(communicator(), path.c_str(), FILE_MODE_WRITE, &fd));
        err = err && check_error(MatView(raw_type(), fd));

        check_error(PetscViewerDestroy(&fd));
        return err;
    }

    bool PetscMatrix::write_matlab(const std::string &path) const {
        PetscViewer fd;

        bool err = check_error(PetscViewerASCIIOpen(communicator(), path.c_str(), &fd));
        err = err && check_error(PetscViewerPushFormat(fd, PETSC_VIEWER_ASCII_MATLAB));
        err = err && check_error(MatView(raw_type(), fd));

        check_error(PetscViewerDestroy(&fd));
        return err;
    }

    void PetscMatrix::copy_from(Mat mat) {
        destroy();
        check_error(MatDuplicate(mat, MAT_COPY_VALUES, &raw_type()));
    }

    void PetscMatrix::copy_to(Mat mat) const { check_error(MatCopy(raw_type(), mat, DIFFERENT_NONZERO_PATTERN)); }

    void PetscMatrix::copy_to(Mat *mat) const { check_error(MatDuplicate(raw_type(), MAT_COPY_VALUES, mat)); }

    void PetscMatrix::transpose() {
        check_error(MatTranspose(raw_type(), MAT_INPLACE_MATRIX, &raw_type()));
        assert(valid());
    }

    void PetscMatrix::transpose(PetscMatrix &result) const {
        if (raw_type() == result.raw_type()) {
            auto s = size();

            if (s.get(0) == s.get(1)) {
                result.transpose();
            } else {
                PetscMatrix temp;
                temp.destroy();

                check_error(MatTranspose(raw_type(), MAT_INITIAL_MATRIX, &temp.raw_type()));
                result.construct(std::move(temp));
            }

            assert(result.valid());
            return;
        }

        result.destroy();
        check_error(MatTranspose(raw_type(), MAT_INITIAL_MATRIX, &result.raw_type()));
        assert(result.valid());
    }

    void PetscMatrix::clear() {
        MPI_Comm comm = communicator();
        destroy();
        MatCreate(comm, &raw_type());
    }

    bool PetscMatrix::is_sparse() const {
        const std::string type_str(type());
        const size_t start = type_str.size() - 5;
        return !(type_str.substr(start, 5) == "dense");
    }

    void PetscMatrix::select(const PetscIndexSet &row_index,
                             const PetscIndexSet &col_index,
                             PetscMatrix &result) const {
        if (col_index.empty()) {
            PetscInt n_rows, n_cols;
            MatGetSize(raw_type(), &n_rows, &n_cols);
            PetscIndexSet all_index(n_cols);
            for (PetscInt i = 0; i < n_cols; ++i) {
                all_index[i] = i;
            }

            select_aux(row_index, all_index, result);

        } else {
            select_aux(row_index, col_index, result);
        }
    }

    void PetscMatrix::select_aux(const PetscIndexSet &row_index,
                                 const PetscIndexSet &col_index,
                                 PetscMatrix &result) const {
        // Mat r = raw_type();
        MPI_Comm comm = communicator();

        SizeType min_col = col_index[0], max_col = col_index[1];
        for (SizeType i : col_index) {
            min_col = std::min(min_col, i);
            max_col = std::max(max_col, i);
        }

        static_assert(std::is_signed<PetscInt>::value, "petsc int must be signed for this to work");
        PetscInt vals[2] = {min_col, -max_col};
        // MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_INT, MPI_MIN, comm);

        this->comm().min(2, vals);
        min_col = vals[0];
        max_col = -vals[1];

        PetscInt global_cols = max_col - min_col + 1;
        PetscInt local_cols = PETSC_DECIDE;
        PetscSplitOwnership(comm, &local_cols, &global_cols);

        unsigned long offsets_in = row_index.size();
        unsigned long offset_out = 0;

        MPI_Exscan(&offsets_in, &offset_out, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);

        offset_out += min_col;

        par_assign_from_local_is(row_index, col_index, min_col, Range(offset_out, offset_out + local_cols), result);
    }

    void PetscMatrix::par_assign_from_local_is(const std::vector<PetscInt> &remote_rows,
                                               const std::vector<PetscInt> &remote_cols,
                                               const PetscInt global_col_offset,
                                               const Range &local_col_range,
                                               PetscMatrix &result) const {
        MPI_Comm comm = communicator();
        Mat &l = result.raw_type();
        const Mat r = raw_type();

        int size;
        MPI_Comm_size(comm, &size);

        int rank;
        MPI_Comm_rank(comm, &rank);

        std::stringstream ss;

        for (auto r : remote_rows) {
            ss << r;
        }

        IS isrow;
        PetscErrorCode ierr;
        UTOPIA_UNUSED(ierr);
        ierr = ISCreateGeneral(comm, remote_rows.size(), &remote_rows[0], PETSC_USE_POINTER, &isrow);

        IS iscol;
        ierr = ISCreateGeneral(comm, remote_cols.size(), &remote_cols[0], PETSC_USE_POINTER, &iscol);

        // TODO(zulianp): maybe make it with collective comms for finiding out if there are off proc entries
        bool has_off_proc_entries = size > 1;

        if (has_off_proc_entries) {
            Mat *l_ptr;

#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 8, 0)
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
            const PetscInt *cols;
            const Scalar *values;

            for (PetscInt row = 0; row < n_rows; ++row) {
                MatGetRow(*l_ptr, row, &n_values, &cols, &values);

                for (PetscInt i = 0; i < n_values; ++i) {
                    MatSetValue(l, rbegin + row, cols[i] - global_col_offset, values[i], INSERT_VALUES);
                }

                MatRestoreRow(*l_ptr, row, &n_values, &cols, &values);
            }

            MatAssemblyBegin(l, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(l, MAT_FINAL_ASSEMBLY);

            MatDestroy(l_ptr);

        } else {
            result.destroy();
#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 8, 0)
            ierr = MatGetSubMatrix(r, isrow, iscol, MAT_INITIAL_MATRIX, &l);
#else
            ierr = MatCreateSubMatrix(r, isrow, iscol, MAT_INITIAL_MATRIX, &l);
#endif  // UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)
        }

        ISDestroy(&isrow);
        ISDestroy(&iscol);
    }

    void PetscMatrix::local_select(const Range &local_row_range,
                                   const Range &local_col_range,
                                   const Range &global_col_range,
                                   PetscMatrix &result) const {
        // PetscErrorCode ierr = 0;

        // Mat &l = result.raw_type();
        // const Mat r = raw_type();

        std::vector<PetscInt> remote_rows;
        remote_rows.reserve(local_row_range.extent());
        for (PetscInt l_row = local_row_range.begin(); l_row < local_row_range.end(); ++l_row) {
            remote_rows.push_back(l_row);
        }

        std::vector<PetscInt> remote_cols;
        remote_cols.reserve(global_col_range.extent());
        for (PetscInt l_col = global_col_range.begin(); l_col < global_col_range.end(); ++l_col) {
            remote_cols.push_back(l_col);
        }

        par_assign_from_local_is(remote_rows, remote_cols, global_col_range.begin(), local_col_range, result);
    }

    void PetscMatrix::select(const Range &global_row_range, const Range &global_col_range, PetscMatrix &result) const {
        const Mat r = raw_type();
        PetscInt global_rows = global_row_range.extent();
        PetscInt local_rows = PETSC_DECIDE;

        PetscInt global_cols = global_col_range.extent();
        PetscInt local_cols = PETSC_DECIDE;

        MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(r));
        PetscSplitOwnership(comm, &local_rows, &global_rows);
        PetscSplitOwnership(comm, &local_cols, &global_cols);

        unsigned long offsets_in[2] = {static_cast<unsigned long>(local_rows), static_cast<unsigned long>(local_cols)};

        unsigned long offset_out[2] = {static_cast<unsigned long>(0), static_cast<unsigned long>(0)};

        MPI_Exscan(offsets_in, offset_out, 2, MPI_UNSIGNED_LONG, MPI_SUM, comm);

        offset_out[0] += global_row_range.begin();
        offset_out[1] += global_col_range.begin();

        local_select(Range(offset_out[0], offset_out[0] + local_rows),
                     Range(offset_out[1], offset_out[1] + local_cols),
                     global_col_range,
                     result);
    }

    PetscMatrix::Scalar PetscMatrix::sum() const {
        Scalar result = static_cast<Scalar>(0);
        MPI_Comm comm = communicator();

        int size = 0;
        MPI_Comm_size(comm, &size);

        if (!is_block() && (size == 1 || !is_mpi())) {
            Vec v;
            MatCreateVecs(raw_type(), nullptr, &v);
            MatGetRowSum(raw_type(), v);
            check_error(VecSum(v, &result));

            VecDestroy(&v);
        } else {
            result = generic_local_reduce(*this, result, utopia::Plus());
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
        }

        return result;
    }

    PetscMatrix::Scalar PetscMatrix::max() const {
        Scalar result = -std::numeric_limits<Scalar>::max();
        MPI_Comm comm = communicator();

        int size = 0;
        MPI_Comm_size(comm, &size);

        if (!is_block() && (size == 1 || !is_mpi())) {
            Vec v;
            MatCreateVecs(raw_type(), nullptr, &v);
            MatGetRowMax(raw_type(), v, nullptr);
            check_error(VecMax(v, nullptr, &result));

            VecDestroy(&v);
        } else {
            result = generic_local_reduce(*this, result, utopia::Max());
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX, comm);
        }

        return result;
    }

    PetscMatrix::Scalar PetscMatrix::min() const {
        Scalar result = std::numeric_limits<Scalar>::max();
        MPI_Comm comm = communicator();

        int size = 0;
        MPI_Comm_size(comm, &size);

        if (!is_block() && (size == 1 || !is_mpi())) {
            Vec v;
            MatCreateVecs(raw_type(), nullptr, &v);
            MatGetRowMin(raw_type(), v, nullptr);
            check_error(VecMin(v, nullptr, &result));

            VecDestroy(&v);
        } else {
            result = generic_local_reduce(*this, result, utopia::Min());
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MIN, comm);
        }

        return result;
    }

    bool PetscMatrix::is_mpi() const {
        static const std::string seq = "seq";
        const std::string str = type();
        return std::search(begin(str), end(str), begin(seq), end(seq)) == str.end();
    }

    void PetscMatrix::build_diag(PetscVector &result) const {
        auto gs = size();

        const bool r_less_than_c = gs.get(0) < gs.get(1);
        const SizeType n = r_less_than_c ? gs.get(0) : gs.get(1);

        if (!result.empty() && n == result.size()) {
            MatGetDiagonal(raw_type(), result.raw_type());
            return;
        }

        result.destroy();

        if (r_less_than_c) {
            MatCreateVecs(raw_type(), nullptr, &result.raw_type());
        } else {
            MatCreateVecs(raw_type(), &result.raw_type(), nullptr);
        }

        check_error(MatGetDiagonal(raw_type(), result.raw_type()));
        result.set_initialized(true);
    }

    void PetscMatrix::col(const SizeType id, PetscVector &result) const {
        auto gs = size();

        result.destroy();
        MatCreateVecs(raw_type(), nullptr, &result.raw_type());

        check_error(MatGetColumnVector(raw_type(), result.raw_type(), id));
        result.set_initialized(true);
    }

    void PetscMatrix::dense_init_diag(MatType dense_type, const PetscVector &diag) {
        MPI_Comm comm = diag.communicator();
        PetscInt local_size = diag.local_size();
        PetscInt global_size = diag.size();

        destroy();

        dense_init(comm, dense_type, local_size, local_size, global_size, global_size);

        check_error(MatZeroEntries(raw_type()));
        check_error(MatDiagonalSet(raw_type(), diag.raw_type(), INSERT_VALUES));
    }

    void PetscMatrix::matij_init_diag(const PetscVector &diag) {
        auto local_size = diag.local_size();
        auto global_size = diag.size();

        destroy();

        matij_init(diag.communicator(), type_override(), local_size, local_size, global_size, global_size, 1, 0);

        check_error(MatZeroEntries(raw_type()));
        check_error(MatDiagonalSet(raw_type(), diag.raw_type(), INSERT_VALUES));
    }

    void PetscMatrix::nest(MPI_Comm comm,
                           PetscInt nr,
                           const IS is_row[],
                           PetscInt nc,
                           const IS is_col[],
                           const Mat a[],
                           const bool use_mat_nest_type) {
        destroy();

        if (use_mat_nest_type) {
            check_error(MatCreateNest(comm, nr, is_row, nc, is_col, a, &raw_type()));
        } else {
            Mat temp = nullptr;

            check_error(MatCreateNest(comm, nr, is_row, nc, is_col, a, &temp));
            check_error(MatConvert(temp, type_override(), MAT_INITIAL_MATRIX, &raw_type()));
            check_error(MatDestroy(&temp));
        }
    }

    void PetscMatrix::diag(const PetscVector &other) {
        // if(is_sparse()) {
        matij_init_diag(other);
        // } else {
        //     dense_init_diag(MATDENSE, other);
        // }
    }

    void PetscMatrix::build_diag(PetscMatrix &result) const {
        MatType type = this->type();
        MPI_Comm comm = communicator();

        result.destroy();
        check_error(MatCreate(comm, &result.raw_type()));
        check_error(MatSetType(result.raw_type(), type));

        const Size gs = size();
        const Size ls = local_size();
        // const bool is_row = gs.get(0) < gs.get(1);

        PetscVector vec;
        build_diag(vec);

        const PetscInt local_size = vec.local_size();
        const PetscInt global_size = vec.size();

        check_error(MatSetSizes(result.raw_type(), local_size, local_size, global_size, global_size));

        // in case it is a sparse format.

        // FIXME handle other cases
        check_error(MatSeqAIJSetPreallocation(result.raw_type(), 1, PETSC_NULL));
        check_error(MatMPIAIJSetPreallocation(result.raw_type(), 1, PETSC_NULL, 0, PETSC_NULL));
        check_error(MatSetUp(result.raw_type()));
        check_error(MatZeroEntries(result.raw_type()));

        check_error(MatDiagonalSet(result.raw_type(), vec.raw_type(), INSERT_VALUES));

        UTOPIA_REPORT_ALLOC("PetscMatrix::build_diag");
    }

    bool PetscMatrix::empty() const {
        Size s = size();
        return s.get(0) <= 0;
    }

    bool PetscMatrix::is_initialized_as(MPI_Comm /*comm*/,
                                        MatType /*dense_type*/,
                                        PetscInt local_rows,
                                        PetscInt local_cols,
                                        PetscInt global_rows,
                                        PetscInt global_cols) {
        if (empty()) {
            return false;
        }

        // TODO(zulianp): : check type and comm

        PetscBool initialized;
        MatAssembled(raw_type(), &initialized);

        if ((initialized != 0u) && (local_rows > 0 && global_cols > 0)) {
            PetscInt m, n;
            MatGetLocalSize(raw_type(), &m, &n);
            initialized = (m == local_rows && n == local_cols) ? PETSC_TRUE : PETSC_FALSE;
        }

        if (initialized != 0u) {
            PetscInt m, n;
            MatGetSize(raw_type(), &m, &n);
            initialized = (m == global_rows && n == global_cols) ? PETSC_TRUE : PETSC_FALSE;
        }

        return initialized == PETSC_TRUE;
    }

    void PetscMatrix::dense_init_values(MPI_Comm comm,
                                        MatType dense_type,
                                        PetscInt local_rows,
                                        PetscInt local_cols,
                                        PetscInt global_rows,
                                        PetscInt global_cols,
                                        Scalar value) {
        if (!is_initialized_as(comm, dense_type, local_rows, local_cols, global_rows, global_cols)) {
            dense_init(comm, dense_type, local_rows, local_cols, global_rows, global_cols);
        }

        const auto r = row_range();
        const PetscInt r_begin = r.begin();
        const PetscInt r_end = r.end();

        const PetscInt computed_global_cols = (global_cols <= 0) ? size().get(1) : global_cols;

        for (PetscInt i = r_begin; i < r_end; ++i) {
            for (PetscInt j = 0; j < computed_global_cols; ++j) {
                MatSetValue(raw_type(), i, j, value, INSERT_VALUES);
            }
        }

        MatAssemblyBegin(raw_type(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(raw_type(), MAT_FINAL_ASSEMBLY);
    }

    void PetscMatrix::dense_init_identity(MPI_Comm comm,
                                          MatType dense_type,
                                          PetscInt local_rows,
                                          PetscInt local_cols,
                                          PetscInt global_rows,
                                          PetscInt global_cols,
                                          Scalar scale_factor) {
        if (!is_initialized_as(comm, dense_type, local_rows, local_cols, global_rows, global_cols)) {
            dense_init(comm, dense_type, local_rows, local_cols, global_rows, global_cols);
        }

        check_error(MatZeroEntries(raw_type()));

        write_lock(utopia::LOCAL);

        const auto r = row_range();
        const PetscInt r_begin = r.begin();

        // otherwise global_cols gives -1, as it should be determined...
        MatGetSize(raw_type(), &global_rows, &global_cols);
        const PetscInt r_end = PetscMin(r.end(), global_cols);

        for (PetscInt i = r_begin; i < r_end; ++i) {
            set(i, i, scale_factor);
        }

        write_unlock(utopia::LOCAL);
    }

    // void PetscMatrix::matij_init_identity(
    //   MPI_Comm comm,
    //   PetscInt local_rows,
    //   PetscInt local_cols,
    //   PetscInt global_rows,
    //   PetscInt global_cols,
    //   Scalar scale_factor)
    // {
    //     matij_init_identity(comm, MATAIJ, local_rows, local_cols, global_rows, global_cols, scale_factor);
    // }

    void PetscMatrix::matij_init_identity(MPI_Comm comm,
                                          MatType sparse_type,
                                          PetscInt local_rows,
                                          PetscInt local_cols,
                                          PetscInt global_rows,
                                          PetscInt global_cols,
                                          Scalar scale_factor) {
        if (!is_initialized_as(comm, sparse_type, local_rows, local_cols, global_rows, global_cols)) {
            matij_init(comm, sparse_type, local_rows, local_cols, global_rows, global_cols, 1, 0);
        }

        MatZeroEntries(raw_type());

        write_lock(utopia::LOCAL);

        const auto r = row_range();
        const auto gs = size();

        const PetscInt r_begin = r.begin();
        const PetscInt r_end = PetscMin(r.end(), gs.get(1));

        for (PetscInt i = r_begin; i < r_end; ++i) {
            set(i, i, scale_factor);
        }

        write_unlock(utopia::LOCAL);

        // MatShift(raw_type(), scale_factor);
    }

    void PetscMatrix::identity(const Scalar &diag) {
        set(0.0);
        shift_diag(diag);
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

    void PetscMatrix::crs(const MatrixLayout &layout,
                          const IndexArray &row_ptr,
                          const IndexArray &col_idx,
                          const ScalarArray &values) {
        // TODO
        // MatCreateSeqAIJWithArrays
        // MatCreateMPIAIJWithArrays
        // MatUpdateMPIAIJWithArrays

        // Se also
        // MatCreateMPIAIJWithSplitArrays
        // MatMPIAIJSetPreallocationCSR

        SizeType n_local_rows = row_ptr.size() - 1;
        assert(n_local_rows == layout.local_size(0));

        assert(layout.comm().size() == 1);  // IMPLEMENT ME for parallel

        IndexArray d_nnz(n_local_rows, 0), o_nnz(n_local_rows, 0);

        for (SizeType r = 0; r < n_local_rows; ++r) {
            d_nnz[r] = row_ptr[r + 1] - row_ptr[r];
        }

        sparse(layout, d_nnz, o_nnz);

        SizeType r_begin = row_range().begin();

        write_lock(utopia::LOCAL);

        for (SizeType r = 0; r < n_local_rows; ++r) {
            for (SizeType idx = row_ptr[r]; idx < row_ptr[r + 1]; ++idx) {
                this->set(r + r_begin, col_idx[idx], values[idx]);
            }
        }

        write_unlock(utopia::LOCAL);
    }

    void PetscMatrix::matij_init(MPI_Comm comm,
                                 MatType type,
                                 PetscInt rows_local,
                                 PetscInt cols_local,
                                 PetscInt rows_global,
                                 PetscInt cols_global,
                                 PetscInt d_nnz,
                                 PetscInt o_nnz) {
        destroy();

        check_error(MatCreate(comm, &raw_type()));
        check_error(MatSetSizes(raw_type(), rows_local, cols_local, rows_global, cols_global));

        check_error(MatSetType(raw_type(), type));

        // MatSetLocalToGlobalMapping HERE (row map and col map)

        check_error(MatSeqAIJSetPreallocation(raw_type(), PetscMax(d_nnz, 1), PETSC_NULL));
        check_error(
            MatMPIAIJSetPreallocation(raw_type(), PetscMax(d_nnz, 1), PETSC_NULL, PetscMax(o_nnz, 1), PETSC_NULL));

        // check_error(
        //     MatISSetPreallocation(raw_type(), PetscMax(d_nnz, 1), PETSC_NULL, PetscMax(o_nnz, 1), PETSC_NULL));

        check_error(MatSetOption(raw_type(), MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE));
        check_error(MatSetOption(raw_type(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE));
        check_error(MatSetOption(raw_type(), MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE));

        check_error(MatZeroEntries(raw_type()));  // TODO not necessary!

        UTOPIA_REPORT_ALLOC("PetscMatrix::matij_init");
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
                                 const IndexArray &d_nnz,
                                 const IndexArray &o_nnz) {
        destroy();

        check_error(MatCreate(comm, &raw_type()));
        check_error(MatSetSizes(raw_type(), rows_local, cols_local, rows_global, cols_global));

        check_error(MatSetType(raw_type(), type));
        check_error(MatSeqAIJSetPreallocation(raw_type(), PETSC_DEFAULT, &d_nnz[0]));
        check_error(MatMPIAIJSetPreallocation(raw_type(), PETSC_DEFAULT, &d_nnz[0], PETSC_DEFAULT, &o_nnz[0]));

        check_error(MatSetOption(raw_type(), MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE));
        check_error(MatSetOption(raw_type(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE));
        check_error(MatSetOption(raw_type(), MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE));

        check_error(MatZeroEntries(raw_type()));

        UTOPIA_REPORT_ALLOC("PetscMatrix::matij_init");
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
                                    PetscInt block_size) {
        destroy();

        check_error(MatCreate(comm, &raw_type()));
        check_error(MatSetSizes(raw_type(), rows_local, cols_local, rows_global, cols_global));
        check_error(MatSetBlockSize(raw_type(), block_size));

        check_error(MatSetType(raw_type(), MATBAIJ));
        check_error(MatSeqBAIJSetPreallocation(raw_type(), block_size, PetscMax(d_nnz, 1), PETSC_NULL));
        check_error(MatMPIBAIJSetPreallocation(
            raw_type(), block_size, PetscMax(d_nnz, 1), PETSC_NULL, PetscMax(o_nnz, 1), PETSC_NULL));

        check_error(MatSetOption(raw_type(), MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE));
        check_error(MatSetOption(raw_type(), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE));
        check_error(MatSetOption(raw_type(), MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE));

        check_error(MatZeroEntries(raw_type()));

        UTOPIA_REPORT_ALLOC("PetscMatrix::mat_baij_init");
    }

    void PetscMatrix::mat_baij_init(MPI_Comm comm,
                                    SizeType rows_local,
                                    SizeType cols_local,
                                    SizeType rows_global,
                                    SizeType cols_global,
                                    const IndexArray &d_nnz,
                                    const IndexArray &o_nnz,
                                    SizeType block_size) {
        destroy();
        check_error(MatCreateBAIJ(comm,
                                  block_size,
                                  rows_local,
                                  cols_local,
                                  rows_global,
                                  cols_global,
                                  -1,
                                  &d_nnz[0],
                                  -1,
                                  &o_nnz[0],
                                  &raw_type()));

        check_error(MatZeroEntries(raw_type()));

        update_mirror();
    }

    bool PetscMatrix::has_nan_or_inf() const {
        int has_nan = 0;
        const Scalar *values;
        const PetscInt *cols;

        PetscInt r_begin, r_end;
        PetscInt n_values = 0;

        PetscInt local_r, local_c;
        MatGetLocalSize(raw_type(), &local_r, &local_c);
        MatGetOwnershipRange(raw_type(), &r_begin, &r_end);

        for (PetscInt row = r_begin; row < r_end; ++row) {
            MatGetRow(raw_type(), row, &n_values, &cols, &values);

            for (PetscInt i = 0; i < n_values; ++i) {
                has_nan = (PetscIsInfOrNanScalar(values[i]) != PetscErrorCode(0) ? 1 : 0);
                if (has_nan != 0) {
                    break;
                }
            }

            MatRestoreRow(raw_type(), row, &n_values, &cols, &values);
            if (has_nan != 0) {
                break;
            }
        }

        if (is_mpi()) {
            MPI_Allreduce(MPI_IN_PLACE, &has_nan, 1, MPI_INT, MPI_MAX, communicator());
        }

        return has_nan > 0;
    }

    void PetscMatrix::inverse(PetscMatrix &result) const {
        assert(!is_sparse());

        PetscMatrix I, L = *this;

        Size gs = size();
        Size ls = local_size();

        I.dense_init(communicator(), type(), ls.get(0), ls.get(1), gs.get(0), gs.get(1));

        check_error(MatZeroEntries(I.raw_type()));

        // initialization
        I.write_lock(utopia::AUTO);
        I.write_unlock(utopia::AUTO);

        check_error(MatShift(I.raw_type(), 1.));

        result.dense_init(communicator(), type(), ls.get(0), ls.get(1), gs.get(0), gs.get(1));

        check_error(MatZeroEntries(result.raw_type()));

        IS isr, isc;
        MatFactorInfo info;

        check_error(MatGetOrdering(L.raw_type(), MATORDERINGNATURAL, &isr, &isc));
        check_error(MatLUFactor(L.raw_type(), isr, isc, &info));
        check_error(MatMatSolve(L.raw_type(), I.raw_type(), result.raw_type()));

        check_error(ISDestroy(&isr));
        check_error(ISDestroy(&isc));
    }

    void PetscMatrix::row_sum(PetscVector &col) const {
        MPI_Comm comm = communicator();

        if (col.is_null() || col.size() != size().get(0)) {
            col.destroy();
            MatCreateVecs(raw_type(), nullptr, &col.raw_type());
        }

        int size = 0;
        MPI_Comm_size(comm, &size);

        if (size == 1 || !is_mpi()) {
            MatGetRowSum(raw_type(), col.raw_type());
            col.set_initialized(true);
        } else {
            reduce_rows(col, *this, 0., utopia::Plus());
        }
    }

    void PetscMatrix::row_max(PetscVector &col) const {
        MPI_Comm comm = communicator();

        if (col.is_null() || col.size() != size().get(0)) {
            col.destroy();
            MatCreateVecs(raw_type(), nullptr, &col.raw_type());
        }

        int size = 0;
        MPI_Comm_size(comm, &size);

        if (size == 1 || !is_mpi()) {
            MatGetRowMax(raw_type(), col.raw_type(), nullptr);
        } else {
            reduce_rows(col, *this, -std::numeric_limits<Scalar>::max(), utopia::Max());
        }
    }

    void PetscMatrix::row_abs_max(PetscVector &col) const {
        MPI_Comm comm = communicator();

        if (col.is_null() || col.size() != size().get(0)) {
            col.destroy();
            MatCreateVecs(raw_type(), nullptr, &col.raw_type());
        }

        int size = 0;
        MPI_Comm_size(comm, &size);

        if (size == 1 || !is_mpi()) {
            MatGetRowMaxAbs(raw_type(), col.raw_type(), nullptr);
        } else {
            reduce_rows(col, *this, 0.0, utopia::AbsMax());
        }
    }

    void PetscMatrix::row_min(PetscVector &col) const {
        MPI_Comm comm = communicator();

        if (col.is_null() || col.size() != size().get(0)) {
            col.destroy();
            MatCreateVecs(raw_type(), nullptr, &col.raw_type());
        }

        int size = 0;
        MPI_Comm_size(comm, &size);

        if (size == 1 || !is_mpi()) {
            MatGetRowMin(raw_type(), col.raw_type(), nullptr);
        } else {
            reduce_rows(col, *this, std::numeric_limits<Scalar>::max(), utopia::Min());
        }
    }

    void PetscMatrix::col_sum(PetscVector &col) const {
        PetscVector temp;
        temp.values(communicator(), col.type(), local_size().get(0), size().get(0), 1.);
        this->transpose_multiply(temp, col);
    }

    void PetscMatrix::multiply(const PetscVector &vec, PetscVector &result) const {
        // handle alias
        if (result.is_alias(vec)) {
            PetscVector x;
            x.copy(vec);
            multiply(x, result);
            return;
        }

        assert(vec.is_consistent());
        assert(result.is_consistent());

        MPI_Comm comm = vec.communicator();

        if (result.is_null()) {
            // MatCreateVecs(raw_type(), nullptr, &result.raw_type());
            create_vecs(nullptr, &result.raw_type());
        } else if (comm != result.communicator()) {
            result.destroy();
            // MatCreateVecs(raw_type(), nullptr, &result.raw_type());
            create_vecs(nullptr, &result.raw_type());
        } else {
            Size gs = size();
            if (gs.get(0) != result.size()) {
                result.destroy();
                // MatCreateVecs(raw_type(), nullptr, &result.raw_type());
                create_vecs(nullptr, &result.raw_type());
            }

            assert(local_size().get(0) == result.local_size());
        }

        check_error(MatMult(raw_type(), vec.raw_type(), result.raw_type()));

        result.update_mirror();
        result.set_initialized(true);

        assert(result.raw_type() != nullptr);
        assert(result.is_consistent());
        // assert(result.same_type(vec)); //FIXME
    }

    void PetscMatrix::transpose_multiply(const PetscVector &vec, PetscVector &result) const {
        if (result.is_alias(vec)) {
            PetscVector temp;
            transpose_multiply(vec, temp);
            result = std::move(temp);
            return;
        }

        MPI_Comm comm = vec.communicator();

        if (result.is_null()) {
            MatCreateVecs(raw_type(), &result.raw_type(), nullptr);
        } else if (comm != result.communicator()) {
            result.destroy();
            MatCreateVecs(raw_type(), &result.raw_type(), nullptr);
        } else {
            Size gs = size();
            Size ls = local_size();
            VecSetSizes(result.raw_type(), ls.get(1), gs.get(1));
        }

        check_error(MatMultTranspose(raw_type(), vec.raw_type(), result.raw_type()));
        assert(result.raw_type() != nullptr);

        result.update_mirror();
        result.set_initialized(true);
    }

    void PetscMatrix::multiply_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const {
        if (v1.raw_type() == result.raw_type() || v2.raw_type() == result.raw_type()) {
            PetscVector temp;
            temp.repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size());
            MatMultAdd(raw_type(), v1.raw_type(), v2.raw_type(), temp.raw_type());
            result = std::move(temp);
        } else {
            result.repurpose(v2.communicator(), v2.type(), v2.local_size(), v2.size());
            MatMultAdd(raw_type(), v1.raw_type(), v2.raw_type(), result.raw_type());
        }

        result.update_mirror();
    }

    void PetscMatrix::transpose_multiply_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const {
        if (v1.raw_type() == result.raw_type() || v2.raw_type() == result.raw_type()) {
            PetscVector temp;
            temp.repurpose(v1.communicator(), v1.type(), v1.local_size(), v1.size());
            MatMultTransposeAdd(raw_type(), v1.raw_type(), v2.raw_type(), temp.raw_type());
            result = std::move(temp);
        } else {
            result.repurpose(v2.communicator(), v2.type(), v2.local_size(), v2.size());
            MatMultTransposeAdd(raw_type(), v1.raw_type(), v2.raw_type(), result.raw_type());
        }

        result.update_mirror();
    }

    /// y = alpha * A * x
    void PetscMatrix::multiply(const Scalar &alpha, const PetscVector &x, PetscVector &y) const {
        multiply(x, y);
        y.scale(alpha);
    }

    /// y = alpha * A^T * x
    void PetscMatrix::transpose_multiply(const Scalar &alpha, const PetscVector &x, PetscVector &y) const {
        transpose_multiply(x, y);
        y.scale(alpha);
    }

    /// y := alpha * A * x + beta * y
    void PetscMatrix::multiply_add(const Scalar &alpha,
                                   const PetscVector &x,
                                   const Scalar &beta,
                                   PetscVector &y) const {
        if (alpha == 1.0 && beta == 1.0) {
            multiply_add(x, y, y);
            return;
        }

        assert(false && "TODO");
    }

    /// y := alpha * A' * x + beta * y
    void PetscMatrix::transpose_multiply_add(const Scalar &alpha,
                                             const PetscVector &x,
                                             const Scalar &beta,
                                             PetscVector &y) const {
        if (alpha == 1.0 && beta == 1.0) {
            transpose_multiply_add(x, y, y);
            return;
        }

        assert(false && "TODO");
    }

    /// y := alpha * op(A) * x + beta * y
    void PetscMatrix::gemv(const bool transpose_A,
                           const Scalar &alpha,
                           const PetscVector &x,
                           const Scalar &beta,
                           PetscVector &y) const {
        if (alpha == 1.0 && beta == 1.0) {
            if (transpose_A) {
                transpose_multiply(x, y);
            } else {
                multiply(x, y);
            }

            return;
        }
        if (beta == 0.0) {
            if (alpha == 0.0) {
                y.set(0.0);
            } else if (transpose_A) {
                transpose_multiply(alpha, x, y);
            } else {
                multiply(alpha, x, y);
            }

            return;
        }

        PetscVector temp;
        temp.copy(y);
        temp.scale(beta);

        multiply(alpha, x, y);
        y.axpy(-1.0, temp);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void PetscMatrix::multiply(const PetscMatrix &mat, PetscMatrix &result) const {
        PetscBool flg;
        // this is very unefficient hack, but still better than fail...
        PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(mat.raw_type()), &flg, MATMPIDENSE, NULL);
        if (flg != 0u) {
            if (mat.raw_type() != result.raw_type() && raw_type() != result.raw_type()) {
                result.destroy();
                Mat temp;
                MatConvert(raw_type(), MATMPIAIJ, MAT_INITIAL_MATRIX, &temp);
                MatMatMult(temp, mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type());
                MatDestroy(&temp);
            } else {
                PetscMatrix temp2;
                temp2.destroy();

                Mat temp;
                MatConvert(raw_type(), MATMPIAIJ, MAT_INITIAL_MATRIX, &temp);
                MatMatMult(temp, mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2.raw_type());
                MatDestroy(&temp);
                result = std::move(temp2);
            }
        } else {
            if (mat.raw_type() != result.raw_type() && raw_type() != result.raw_type()) {
                result.destroy();
                MatMatMult(raw_type(), mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type());
            } else {
                PetscMatrix temp;
                temp.destroy();
                MatMatMult(raw_type(), mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.raw_type());
                result = std::move(temp);
            }
        }

        result.update_mirror();
        assert(result.valid());
    }

    void PetscMatrix::transpose_multiply(const PetscMatrix &mat, PetscMatrix &result) const {
        if (mat.raw_type() != result.raw_type() && raw_type() != result.raw_type()) {
            result.destroy();
            check_error(
                MatTransposeMatMult(raw_type(), mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()));
            assert(result.valid());
        } else {
            PetscMatrix temp;
            temp.destroy();
            check_error(
                MatTransposeMatMult(raw_type(), mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.raw_type()));
            result.construct(std::move(temp));
            assert(result.valid());
        }

        result.update_mirror();
    }

    void PetscMatrix::multiply_transpose(const PetscMatrix &mat, PetscMatrix &result) const {
        m_utopia_warning("> FIXME MatMatTransposeMult does not work in parallel, prepare work around");

        if (mat.raw_type() != result.raw_type() && raw_type() != result.raw_type()) {
            result.destroy();
            check_error(
                MatMatTransposeMult(raw_type(), mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()));
            assert(result.valid());
        } else {
            PetscMatrix temp;
            temp.destroy();
            check_error(
                MatMatTransposeMult(raw_type(), mat.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp.raw_type()));
            result = std::move(temp);
            assert(result.valid());
        }

        result.update_mirror();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void PetscMatrix::axpy(const Scalar &a, const PetscMatrix &x) {
        UTOPIA_TRACE_SCOPE("PetscMatrix::axpy");
        check_error(MatAXPY(raw_type(), a, x.raw_type(), DIFFERENT_NONZERO_PATTERN));
    }

    void PetscMatrix::axpy_subset(const Scalar &a, const PetscMatrix &x) {
        UTOPIA_TRACE_SCOPE("PetscMatrix::axpy_subset");
        check_error(MatAXPY(raw_type(), a, x.raw_type(), SUBSET_NONZERO_PATTERN));
    }

    void PetscMatrix::convert_to_mat_baij(const PetscInt block_size, PetscMatrix &output) {
        auto ls = local_size();
        auto gs = size();

        output.mat_baij_init(
            communicator(), ls.get(0), ls.get(1), gs.get(0), gs.get(1), PETSC_DEFAULT, PETSC_DEFAULT, block_size);

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0)
        output.write_lock(utopia::AUTO);
        output.write_unlock(utopia::AUTO);
#endif  // UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0)

        check_error(MatCopy(raw_type(), output.raw_type(), DIFFERENT_NONZERO_PATTERN));
    }

    void PetscMatrix::convert_to_mat_baij(const PetscInt block_size) {
        auto ls = local_size();
        auto gs = size();

        PetscMatrix temp;
        convert_to_mat_baij(block_size, temp);
        *this = std::move(temp);
    }

    bool PetscMatrix::is_assembled() const {
        if (empty()) return false;

        PetscBool assembled = PETSC_FALSE;
        check_error(MatAssembled(raw_type(), &assembled));
        return assembled == PETSC_TRUE;
    }

    // testing MATAIJCUSPARSE,MATSEQAIJCUSPARSE
    bool PetscMatrix::PetscMatrix::is_cuda() const {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(raw_type()), MATAIJCUSPARSE, &match);
        if (match == PETSC_TRUE) {
            return true;
        }

        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(raw_type()), MATSEQAIJCUSPARSE, &match);
        return match == PETSC_TRUE;
    }

    bool PetscMatrix::is_block(Mat mat) {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(mat), MATBAIJ, &match);
        if (match == PETSC_TRUE) {
            return true;
        }

        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(mat), MATSEQBAIJ, &match);
        if (match == PETSC_TRUE) {
            return true;
        }

        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(mat), MATMPIBAIJ, &match);
        return match == PETSC_TRUE;
    }

    bool PetscMatrix::is_block() const { return is_block(raw_type()); }

    VecType PetscMatrix::compatible_cuda_vec_type() const {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(raw_type()), MATAIJCUSPARSE, &match);

        if (match != 0u) {
            return VECMPICUDA;
        }
        return VECSEQCUDA;
    }

    bool PetscMatrix::create_vecs(Vec *x, Vec *y) const {
        MatCreateVecs(raw_type(), x, y);
        return true;
    }

    bool PetscMatrix::has_type(MatType type) const {
        PetscBool match = PETSC_FALSE;
        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(raw_type()), type, &match);
        return match == PETSC_TRUE;
    }

    bool PetscMatrix::same_type(const PetscMatrix &other) const { return has_type(other.type()); }

    void PetscMatrix::set_zero_rows(const PetscIndexSet &idx, const Scalar &diag) {
        PetscBool val = PETSC_TRUE;

        MatGetOption(raw_type(), MAT_KEEP_NONZERO_PATTERN, &val);
        MatSetOption(raw_type(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

        check_error(MatZeroRows(raw_type(), idx.size(), &idx[0], diag, nullptr, nullptr));

        MatSetOption(raw_type(), MAT_KEEP_NONZERO_PATTERN, val);
    }

    bool PetscMatrix::equals(const PetscMatrix &other, const Scalar &tol) const {
        PetscMatrix diff = other;
        diff.axpy(-1.0, *this);
        return diff.norm_infty() <= tol;
    }

    void PetscMatrix::diagonal_block(PetscMatrix &other) const {
        // reference count on the returned matrix is not incremented and it is used as part of the containing MPI
        // Mat's normal operation
        Mat M;
        check_error(MatGetDiagonalBlock(raw_type(), &M));
        other.copy_from(M);
        assert(other.valid());
    }

    void PetscMatrix::diag_scale_right(const PetscVector &diag) {
        check_error(MatDiagonalScale(raw_type(), nullptr, diag.raw_type()));
    }

    void PetscMatrix::diag_scale_left(const PetscVector &diag) {
        check_error(MatDiagonalScale(raw_type(), diag.raw_type(), nullptr));
    }

    void PetscMatrix::convert_from(const Mat &mat) { copy_from(mat); }

    void PetscMatrix::convert_to(Mat &mat) const { copy_to(mat); }

    /// C := alpha * A * B
    void PetscMatrix::multiply(const Scalar &alpha, const PetscMatrix &B, PetscMatrix &C) const {
        multiply(B, C);
        C.scale(alpha);
    }

    /// C := op(A) * op(B)
    void PetscMatrix::multiply(const bool transpose_A,
                               const bool transpose_B,
                               const PetscMatrix &B,
                               PetscMatrix &C) const {
        if (!transpose_A && !transpose_B) {
            multiply(B, C);
            return;
        }

        if (transpose_A && !transpose_B) {
            transpose_multiply(B, C);
            return;
        }

        if (!transpose_A && transpose_B) {
            multiply_transpose(B, C);
            return;
        }

        B.multiply(*this, C);
        C.transpose();
    }

    bool PetscMatrix::valid() const {
        // PetscTruth       flg;
        // MatValid(mat,(PetscTruth*)&flg);
        PetscBool flg;
        MatAssembled(raw_type(), &flg);
        return (flg != 0u) && (type() != nullptr);
    }

    void PetscMatrix::assign(const Range & /*row_range*/, const Range & /*col_range*/, const PetscMatrix & /*block*/) {
        assert(false && "IMPLEMENT ME");
    }

    PetscMatrix::SizeType PetscMatrix::global_nnz() const {
        if (empty()) {
            return 0;
        }

        MatInfo info;
        MatGetInfo(raw_type(), MAT_GLOBAL_SUM, &info);
        return info.nz_used;
    }

    PetscMatrix::SizeType PetscMatrix::local_nnz() const {
        if (empty()) {
            return 0;
        }

        MatInfo info;
        MatGetInfo(raw_type(), MAT_LOCAL, &info);
        return info.nz_used;
    }

    PetscMatrix::SizeType PetscMatrix::nnz(const Scalar tol) const {
        if (empty()) {
            return 0;
        }

        SizeType ret = 0;
        // FIXME reading indices is not necessary
        this->read([tol, &ret](const SizeType &, const SizeType &, const Scalar &v) { ret += PetscAbs(v) > tol; });

        return ret;
    }

    void PetscMatrix::shift_diag(const PetscVector &d) {
        if (empty()) {
            diag(d);
        } else {
            check_error(MatDiagonalSet(raw_type(), d.raw_type(), ADD_VALUES));
        }
    }

    void PetscMatrix::set_diag(const PetscVector &d) {
        check_error(MatDiagonalSet(raw_type(), d.raw_type(), INSERT_VALUES));
    }

    void PetscMatrix::set(const Scalar &val) {
        if (val == 0.0) {
            scale(0.0);
            return;
        }

        if (is_sparse()) {
            assert(false);
            m_utopia_error("trying to set a sparse matrix to a dense value");
            Utopia::Abort();
        } else {
            transform_values([val](const Scalar &) -> Scalar { return val; });
        }
    }

    void PetscMatrix::update_mirror() { comm_.set(communicator()); }

    void PetscMatrix::convert_to_scalar_matrix() {
        if (this->is_block()) {
            PetscMatrix temp;
            temp.destroy();
            MatConvert(this->raw_type(), MATAIJ, MAT_INITIAL_MATRIX, &temp.raw_type());
            *this = std::move(temp);
        }
    }

    void PetscMatrix::convert_to_scalar_matrix(PetscMatrix &scalar_matrix) {
        if (this->is_block()) {
            scalar_matrix.destroy();
            MatConvert(this->raw_type(), MATAIJ, MAT_INITIAL_MATRIX, &scalar_matrix.raw_type());
        } else {
            scalar_matrix = *this;
        }
    }

    // https://petsc.org/release/docs/manualpages/Mat/MatCreateMPIAIJWithSplitArrays/

    void PetscMatrix::wrap(MPI_Comm comm,
                           const PetscInt rows_local,
                           const PetscInt cols_local,
                           const PetscInt rows_global,
                           const PetscInt cols_global,
                           PetscInt *rowptr,
                           PetscInt *colidx,
                           PetscScalar *values,
                           std::function<void()> destroy_callback) {
        destroy();

        int size;
        MPI_Comm_size(comm, &size);

        if (size == 1) {
            check_error(MatCreateSeqAIJWithArrays(comm, rows_global, cols_global, rowptr, colidx, values, &raw_type()));

        } else {
            check_error(MatCreateMPIAIJWithArrays(
                comm, rows_local, cols_local, rows_global, cols_global, rowptr, colidx, values, &raw_type()));
        }

        this->destroy_callback = destroy_callback;
        update_mirror();
    }

    void PetscMatrix::set_block_size(const int block_size) {
        // int prev_block_size;
        // MatGetBlockSize(raw_type(), &prev_block_size);
        // printf("prev_block_size=%d\n", prev_block_size);
        MatSetBlockSize(raw_type(), block_size);
    }

    void PetscMatrix::ghosts(IndexArray &ret) const {
        PetscInt nghosts;
        const PetscInt *ghosts = nullptr;
        MatGetGhosts(this->raw_type(), &nghosts, &ghosts);

        ret.resize(nghosts);
        std::copy(ghosts, ghosts + nghosts, ret.begin());
    }

    void PetscMatrix::lump() {
        PetscCrsView d, o;
        views_host(*this, d, o);

        SizeType rows = d.rows();
        for (SizeType r = 0; r < rows; r++) {
            Scalar row_sum = 0;

            {
                auto row = o.row(r);
                for (SizeType k = 0; k < row.length; k++) {
                    Scalar Aij = row.value(k);
                    row_sum += Aij;
                    row.value(k) = 0;
                }
            }

            {
                auto row = d.row(r);
                SizeType k_diag = -1;
                for (SizeType k = 0; k < row.length; k++) {
                    SizeType c = row.colidx(k);
                    Scalar Aij = row.value(k);

                    row_sum += Aij;
                    row.value(k) = 0;

                    if (c == r) {
                        k_diag = k;
                    }
                }

                assert(k_diag >= 0);

                if (k_diag >= 0) {
                    row.value(k_diag) = row_sum;
                }
            }
        }
    }

}  // namespace utopia

// TODO(zulianp):
// class CompatibleMatPair {
// public:
//     CompatibleMatPair(const MPI_Comm comm, const Mat &left, const Mat &right)
//     {
//         const bool left_is_sparse  = is_sparse(left);
//         const bool right_is_sparse = is_sparse(right);

//         must_destroy_left_  = false;
//         must_destroy_right_ = false;

//         MatType common_type;
//         if(left_is_sparse != right_is_sparse) {
//             if(left_is_sparse) {
//                 MatCreate(comm, &left_);
//                 PetscBackend::check_error( MatGetType(right, &common_type) );
//                 PetscBackend::check_error( MatConvert(left, common_type, MAT_INITIAL_MATRIX, &left_) );
//                 right_ = right;
//                 must_destroy_left_  = true;

//             } else {
//                 MatCreate(comm, &right_);
//                 PetscBackend::check_error( MatGetType(left, &common_type) );
//                 PetscBackend::check_error( MatConvert(right, common_type, MAT_INITIAL_MATRIX, &right_) );
//                 left_ = left;
//                 must_destroy_right_ = true;
//             }
//         } else {
//             left_  = left;
//             right_ = right;
//         }
//     }

//     inline const Mat &left()
//     {
//         return left_;
//     }

//     inline const Mat &right()
//     {
//         return right_;
//     }

//     ~CompatibleMatPair()
//     {
//         if(must_destroy_left_) MatDestroy(&left_);
//         if(must_destroy_right_) MatDestroy(&right_);
//     }

// private:

//     static bool is_sparse(const Mat &mat)
//     {
//         MatType type;
//         MatGetType(mat, &type);
//         const std::string type_str(type);
//         const size_t start = type_str.size() - 5;
//         return !(type_str.substr(start, 5) == "dense");
//     }

//     Mat left_;
//     Mat right_;
//     bool must_destroy_left_;
//     bool must_destroy_right_;
// };
