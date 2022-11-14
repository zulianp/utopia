#ifndef UTOPIA_DISCRETIZATION_MANAGER_HPP
#define UTOPIA_DISCRETIZATION_MANAGER_HPP

#include "utopia_fe_Core.hpp"

#include "utopia_kokkos_FEBase.hpp"

#include "utopia_Field.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_kokkos_Field.hpp"

#include <memory>
#include <string>

namespace utopia {
    namespace kokkos {
        template <class FunctionSpace,
                  class FE,
                  class MatrixAccumulator_ = DefaultView<typename FE::Scalar>,
                  class VectorAccumulator_ = DefaultView<typename FE::Scalar>,
                  class ScalarAccumulator_ = DefaultView<typename FE::Scalar>>
        class Discretization : public Configurable {
        public:
            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using Field = utopia::Field<FunctionSpace>;
            using FEField = utopia::kokkos::Field<FE>;
            using FE_ptr = std::shared_ptr<FE>;
            using FunctionSpace_ptr = std::shared_ptr<FunctionSpace>;

            using MatrixAccumulator = MatrixAccumulator_;
            using VectorAccumulator = VectorAccumulator_;
            using ScalarAccumulator = ScalarAccumulator_;

            class Part {
            public:
                Part(std::string name = "") : name(std::move(name)) {}
                std::string name;
                inline bool is_all() const { return name.empty(); }
            };

            static Part all() {
                static Part part;
                return part;
            }
            ////////////////////////////////////////////////////////////////////////////////////

            virtual void create(std::vector<FE_ptr> &fe, int order, const Part &part = all()) = 0;
            virtual void create_on_boundary(FE &fe, int order, const Part &part = all()) = 0;

            ////////////////////////////////////////////////////////////////////////////////////

            virtual void convert_field(const Field &in,
                                       std::vector<std::shared_ptr<FEField>> &out,
                                       const Part &part = all()) = 0;

            virtual void convert_field(const std::vector<std::shared_ptr<FEField>> &in,
                                       Field &out,
                                       const Part &part = all()) = 0;

            ////////////////////////////////////////////////////////////////////////////////////

            virtual void global_to_local(const Vector &vector,
                                         const std::vector<VectorAccumulator> &element_vectors,
                                         const Part &part = all(),
                                         const int comp = 0) = 0;

            virtual void local_to_global(const std::vector<MatrixAccumulator> &acc,
                                         AssemblyMode mode,
                                         Matrix &mat,
                                         const Part &part = all()) = 0;

            virtual void local_to_global(const std::vector<VectorAccumulator> &acc,
                                         AssemblyMode mode,
                                         Vector &vec,
                                         const Part &part = all()) = 0;

            virtual void local_to_global_on_boundary(const std::vector<VectorAccumulator> &acc,
                                                     AssemblyMode mode,
                                                     Vector &vec,
                                                     const Part &part = all()) = 0;

            ////////////////////////////////////////////////////////////////////////////////
            // Single acc (Retro comp tobe removed)
            virtual void global_to_local(const Vector &vector,
                                         const VectorAccumulator &element_vectors,
                                         const Part &part = all(),
                                         const int comp = 0) {
                global_to_local(vector, {element_vectors}, part, comp);
            }

            virtual void local_to_global(const MatrixAccumulator &acc,
                                         AssemblyMode mode,
                                         Matrix &mat,
                                         const Part &part = all()) {
                local_to_global({acc}, mode, mat, part);
            }

            virtual void local_to_global(const VectorAccumulator &acc,
                                         AssemblyMode mode,
                                         Vector &vec,
                                         const Part &part = all()) {
                local_to_global({acc}, mode, vec, part);
            }

            virtual void local_to_global_on_boundary(const VectorAccumulator &acc,
                                                     AssemblyMode mode,
                                                     Vector &vec,
                                                     const Part &part = all()) {
                local_to_global_on_boundary({acc}, mode, vec, part);
            }

            Discretization(const FunctionSpace_ptr &space) : space_(space) {}
            inline FunctionSpace_ptr space() const { return space_; }

            virtual ~Discretization() = default;

        private:
            FunctionSpace_ptr space_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_DISCRETIZATION_MANAGER_HPP
