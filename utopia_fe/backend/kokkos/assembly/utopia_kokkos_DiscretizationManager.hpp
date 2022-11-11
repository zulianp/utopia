#ifndef UTOPIA_DISCRETIZATION_MANAGER_HPP
#define UTOPIA_DISCRETIZATION_MANAGER_HPP

#include "utopia_fe_Core.hpp"

#include "utopia_kokkos_FEBase.hpp"

#include "utopia_Field.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include <string>

namespace utopia {
    namespace kokkos {
        template <class FunctionSpace,
                  class FE,
                  class MatrixAccumulator = DefaultView<typename FE::Scalar>,
                  class VectorAccumulator = DefaultView<typename FE::Scalar>,
                  class ScalarAccumulator = DefaultView<typename FE::Scalar>>
        class DiscretizationManager : public Configurable {
        public:
            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using Field = utopia::Field<FunctionSpace>;
            using FEField = utopia::kokkos::Field<FE>;

            class Part {
            public:
                Part(std::string name = "") : name(std::move(name)) {}
                std::string name;
            };

            static const Part all;

            void convert_field(const Field &in, FEField &out) = 0;
            void convert_field(const FEField &in, Field &out) = 0;

            void global_to_local(const FunctionSpace &space,
                                 const Vector &vector,
                                 VectorAccumulator &element_vectors,
                                 const int n_comp = 1);

            void create_fe(const FunctionSpace &space, FE &fe, int order, const Part &part = all) = 0;

            void create_fe_on_boundary(const FunctionSpace &space, FE &fe, int order, const Part &part = all) = 0;

            void local_to_global(const FunctionSpace &space,
                                 MatrixAccumulator &acc,
                                 AssemblyMode mode,
                                 Matrix &mat,
                                 const Part &part = all) = 0;

            void local_to_global(const FunctionSpace &space,
                                 VectorAccumulator &acc,
                                 AssemblyMode mode,
                                 Vector &vec,
                                 const Part &part = all) = 0;

            void local_to_global_on_boundary(const FunctionSpace &space,
                                             VectorAccumulator &acc,
                                             AssemblyMode mode,
                                             Vector &vec,
                                             const Part &part = all) = 0;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_DISCRETIZATION_MANAGER_HPP
