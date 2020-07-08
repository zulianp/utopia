#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#ifdef WITH_PETSC

#include "utopia_petsc.hpp"
#include "utopia_petsc_trilinos_ConvertTensor.hpp"
#include "utopia_trilinos.hpp"

// FIXME do not use this!!! It will not be portable in the future
#include "utopia_trilinos_Each_impl.hpp"

#include <numeric>
#include <vector>

namespace utopia {

    void ConvertTensor<TpetraVector, PetscVector, 1, TRILINOS, PETSC>::apply(const TpetraVector &in, PetscVector &out) {
        PetscCommunicator comm(in.comm().raw_comm());
        out.zeros(layout(comm, in.local_size(), in.size()));

        // FIXME do not use this!!! It will not be portable in the future

        auto r = range(in);

        Read<TpetraVector> r_in(in);
        Write<PetscVector> w_out(out);

        for (auto i = r.begin(); i < r.end(); ++i) {
            out.set(i, in.get(i));
        }
    }

    void ConvertTensor<PetscVector, TpetraVector, 1, PETSC, TRILINOS>::apply(const PetscVector &in, TpetraVector &out) {
        TrilinosCommunicator comm(in.comm().raw_comm());
        out.zeros(layout(comm, in.local_size(), in.size()));

        // FIXME do not use this!!! It will not be portable in the future

        auto r = range(in);

        Read<PetscVector> r_in(in);
        Write<TpetraVector> w_out(out);

        for (auto i = r.begin(); i < r.end(); ++i) {
            out.set(i, in.get(i));
        }
    }

    void ConvertTensor<TpetraMatrix, PetscMatrix, 2, TRILINOS, PETSC>::apply(const TpetraMatrix &from,
                                                                             PetscMatrix &to) {
        using Scalar = typename Traits<TpetraMatrix>::Scalar;
        using SizeType = typename Traits<TpetraMatrix>::SizeType;

        // FIXME do not use this!!! It will not be portable in the future

        auto ls = local_size(from);
        auto n_row_local = ls.get(0);
        std::vector<int> nnzxrow(n_row_local, 0);
        auto r = row_range(from);

        TpetraMatrixEach::apply_read(
            from, [&nnzxrow, &r](const SizeType &i, const SizeType &, const Scalar &) { ++nnzxrow[i - r.begin()]; });

        auto nnz = *std::max_element(nnzxrow.begin(), nnzxrow.end());

        // FIXME use nnzxrow instead
        to.sparse(layout(from), nnz, nnz);

        {
            Write<PetscMatrix> w_t(to);
            TpetraMatrixEach::apply_read(
                from, [&to](const SizeType &i, const SizeType &j, const Scalar &val) { to.set(i, j, val); });
        }

        assert(size(from) == size(to));
        assert(local_size(from) == local_size(to));
    }

    void ConvertTensor<PetscMatrix, TpetraMatrix, 2, PETSC, TRILINOS>::apply(const PetscMatrix &from,
                                                                             TpetraMatrix &to) {
        using Scalar = typename Traits<PetscMatrix>::Scalar;
        using SizeType = typename Traits<PetscMatrix>::SizeType;

        // FIXME do not use this!!! It will not be portable in the future

        auto ls = local_size(from);
        auto n_row_local = ls.get(0);
        std::vector<int> nnzxrow(n_row_local, 0);
        auto r = row_range(from);

        from.read([&nnzxrow, &r](const SizeType &i, const SizeType &, const Scalar &) { ++nnzxrow[i - r.begin()]; });

        auto nnz = *std::max_element(nnzxrow.begin(), nnzxrow.end());

        // FIXME use nnzxrow instead
        to.sparse(layout(from), nnz, nnz);

        {
            Write<TpetraMatrix> w_t(to);
            from.read([&to](const SizeType &i, const SizeType &j, const Scalar &val) { to.set(i, j, val); });
        }

        assert(size(from) == size(to));
        assert(local_size(from) == local_size(to));
    }

}  // namespace utopia

#endif
#endif
