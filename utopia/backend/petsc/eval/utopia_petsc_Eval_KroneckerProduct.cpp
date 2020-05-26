#include "utopia_petsc_Eval_KroneckerProduct.hpp"

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_petsc_Communicator.hpp"

#include <mpi.h>
#include <vector>

namespace utopia {

    template <class Matrix, class Vector>
    void EvalKroneckerProduct<Matrix, Vector, PETSC>::apply(const Vector &left, const Vector &right, Matrix &result) {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        MPI_Comm comm = left.comm().get();

        int n_procs = 0, rank = 0;
        MPI_Comm_size(comm, &n_procs);
        MPI_Comm_rank(comm, &rank);

        const SizeType lsize = left.size();
        const SizeType rsize = right.size();

        const Size result_size = {lsize, rsize};

        const Scalar *right_array = nullptr;
        PetscErrorCode err;
        UTOPIA_UNUSED(err);
        err = VecGetArrayRead(right.raw_type(), &right_array);

        const Scalar *left_array = nullptr;
        err = VecGetArrayRead(left.raw_type(), &left_array);

        const Range r_range = range(right);
        const Range l_range = range(left);
        const PetscInt n = rsize;

        std::vector<Scalar> recvbuf(n, 0);

        // not very efficient but good enough for the moment
        int is_evenly_distributed = static_cast<int>(n == r_range.extent() * n_procs);
        MPI_Allreduce(MPI_IN_PLACE, &is_evenly_distributed, 1, MPI_INT, MPI_MIN, comm);

        if (is_evenly_distributed != 0) {
            MPI_Allgather(right_array,
                          r_range.extent(),
                          MPIType<Scalar>::value(),
                          &recvbuf[0],
                          r_range.extent(),
                          MPIType<Scalar>::value(),
                          comm);
        } else {
            const long n_values = r_range.extent();
            std::vector<long> n_values_x_proc(n_procs);
            std::vector<long> offsets(n_procs);
            n_values_x_proc[rank] = n_values;

            MPI_Allgather(&n_values, 1, MPIType<long>::value(), &n_values_x_proc[0], 1, MPIType<long>::value(), comm);

            offsets[0] = 0;
            for (int r = 1; r != n_procs; ++r) {
                offsets[r] = offsets[r - 1] + n_values_x_proc[r - 1];
            }

            std::copy(right_array, right_array + n_values, recvbuf.begin() + offsets[rank]);
            std::vector<MPI_Request> requests((n_procs - 1) * 2);

            SizeType req_index = 0;
            for (int r = 0; r < n_procs; ++r) {
                if (r == rank) {
                    continue;
                }

                MPI_Isend(right_array, r_range.extent(), MPIType<Scalar>::value(), r, r, comm, &requests[req_index++]);
                MPI_Irecv(&recvbuf[offsets[r]],
                          n_values_x_proc[r],
                          MPIType<Scalar>::value(),
                          r,
                          rank,
                          comm,
                          &requests[req_index++]);
            }

            MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
        }

        err = VecRestoreArrayRead(right.raw_type(), &right_array);

        result.dense_init(
            comm, result.type_override(), l_range.extent(), PETSC_DECIDE, result_size.get(0), result_size.get(1));

        result.write_lock(LOCAL);

        for (SizeType i = l_range.begin(); i != l_range.end(); ++i) {
            const Scalar l_value = left_array[i - l_range.begin()];

            for (SizeType j = 0; j != n; ++j) {
                const Scalar r_value = recvbuf.at(j);

                MatSetValue(result.raw_type(), i, j, l_value * r_value, INSERT_VALUES);
            }
        }

        result.write_unlock(LOCAL);

        VecRestoreArrayRead(left.raw_type(), &left_array);
    }

    template class EvalKroneckerProduct<PetscMatrix, PetscVector, PETSC>;
}  // namespace utopia
