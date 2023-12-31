#ifndef UTOPIA_TEMP_HPP
#define UTOPIA_TEMP_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

#include "utopia_Expression.hpp"
#include "utopia_MatChop.hpp"

#include <iostream>
#include <type_traits>

namespace utopia {

    template <class Derived>
    void set_zero_rows(Tensor<Derived, 2> &w,
                       const typename Traits<Derived>::IndexSet &index,
                       const typename Traits<Derived>::Scalar diag) {
        w.derived().set_zero_rows(index, diag);
    }

    template <class Derived>
    void set(Tensor<Derived, 1> &v,
             const typename Traits<Derived>::IndexSet &index,
             const typename Traits<Derived>::Scalar value) {
        auto &vec = v.derived();

        Write<Derived> write_lock(vec);

        for (auto &idx : index) {
            vec.set(idx, value);
        }
    }

    template <class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class ApplyEqualityConstraintsToSystem {
    public:
        template <class IndexSetT>
        void apply(Matrix &, Vector &, Vector &, const IndexSetT &) {
            static_assert(Backend < HOMEMADE, "implement for specific backend");
        }
    };

    // FIXME not sure how to name this one: change name to apply_equality_constraints_to_system
    template <class MatDerived, class VecDerived>
    void apply_BC_to_system(Tensor<MatDerived, 2> &A,
                            Tensor<VecDerived, 1> &x,
                            Tensor<VecDerived, 1> &rhs,
                            const typename Traits<MatDerived>::IndexSet &constrained_idx) {
        ApplyEqualityConstraintsToSystem<MatDerived, VecDerived>::apply(
            A.derived(), x.derived(), rhs.derived(), constrained_idx);
    }

    template <class MatDerived, class VecDerived>
    void apply_equality_constraints_to_system(Tensor<MatDerived, 2> &A,
                                              Tensor<VecDerived, 1> &x,
                                              Tensor<VecDerived, 1> &rhs,
                                              const typename Traits<MatDerived>::IndexSet &constrained_idx) {
        ApplyEqualityConstraintsToSystem<MatDerived, VecDerived>::apply(
            A.derived(), x.derived(), rhs.derived(), constrained_idx);
    }

    template <class Matrix, class Vector>
    void set_zero_rows(Tensor<Matrix, 2> &w, const Tensor<Vector, 1> &indicator, const double diag = 0.) {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        // FIXME maybe use array type once is available
        using IndexSet = typename Traits<Vector>::IndexSet;

        IndexSet index;
        // index.reserve(local_size(indicator).get(0));

        // {
        //     each_read(indicator.derived(), [&index](const SizeType i, const Scalar value) {
        //         if (value == 1.) {
        //             index.push_back(i);
        //         }
        //     });
        // }

        {
            indicator.derived().read([&index](const SizeType i, const Scalar &value) {
                if (value == 1.) {
                    index.push_back(i);
                }
            });
        }

        set_zero_rows(w, index, diag);
    }

    template <class Vector, int Backend = Traits<Vector>::Backend>
    class EvalVecUniqueSortSerial {
    public:
        static void apply(const Tensor<Vector, 1> & /*x*/,
                          Tensor<Vector, 1> & /*sorted */,
                          const int /*used_values */) {
            static_assert(Traits<Vector>::Backend == PETSC,
                          "EvalVecUniqueSortSerial implemented just for petsc backend.");
        }
    };

    template <class Vector>
    void vec_unique_sort_serial(const Tensor<Vector, 1> &x, Tensor<Vector, 1> &sorted, const int used_values = -1) {
        EvalVecUniqueSortSerial<Vector>::apply(x, sorted, used_values);
    }

}  // namespace utopia

#endif  // UTOPIA_TEMP_HPP
