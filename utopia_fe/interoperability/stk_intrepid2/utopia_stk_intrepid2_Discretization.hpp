#ifndef UTOPIA_STK_INTREPID2_DISCRETIZATION_HPP
#define UTOPIA_STK_INTREPID2_DISCRETIZATION_HPP

#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_kokkos_Discretization.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include <memory>

namespace utopia {

    using stk_FS_t = utopia::stk::FunctionSpace;
    using stk_FE_t = utopia::intrepid2::FE<utopia::Traits<utopia::stk::FunctionSpace>::Scalar>;

    template <>
    class Discretization<stk_FS_t, stk_FE_t> final : public utopia::kokkos::Discretization<stk_FS_t, stk_FE_t> {
    public:
        using FunctionSpace = stk_FS_t;
        using FE = stk_FE_t;
        using Super = utopia::kokkos::Discretization<FunctionSpace, FE>;
        using Part = typename Super::Part;

        using VectorAccumulator = typename Super::VectorAccumulator;
        using MatrixAccumulator = typename Super::MatrixAccumulator;
        using ScalarAccumulator = typename Super::ScalarAccumulator;

        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Scalar = typename Traits<FunctionSpace>::Scalar;

        Discretization(const std::shared_ptr<FunctionSpace> &space, const std::shared_ptr<FE> &fe);
        ~Discretization();

        void read(Input &in) override;

        ////////////////////////////////////////////////////////////////////////////////////

        void create(std::vector<FE_ptr> &fe, int order, const Part &part = Super::all()) override;
        void create_on_boundary(std::vector<FE_ptr> &fe, int order, const Part &part = Super::all()) override;

        ////////////////////////////////////////////////////////////////////////////////////

        void convert_field(const Field &in,
                           std::vector<std::shared_ptr<FEField>> &out,
                           const Part &part = Super::all()) override;

        void convert_field(const std::vector<std::shared_ptr<FEField>> &in,
                           Field &out,
                           const Part &part = Super::all()) override;

        ////////////////////////////////////////////////////////////////////////////////////

        void global_to_local(const Vector &vector,
                             std::vector<VectorAccumulator> &element_vectors,
                             const Part &part = Super::all(),
                             const int comp = 0) override;

        // Local to global

        void local_to_global(const std::vector<MatrixAccumulator> &acc,
                             AssemblyMode mode,
                             Matrix &mat,
                             const Part &part = Super::all()) override;

        void local_to_global(const std::vector<VectorAccumulator> &acc,
                             AssemblyMode mode,
                             Vector &vec,
                             const Part &part = Super::all()) override;

        void local_to_global_on_boundary(const std::vector<VectorAccumulator> &acc,
                                         AssemblyMode mode,
                                         Vector &vec,
                                         const Part &part = Super::all()) override;
        class Impl;

    private:
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_DISCRETIZATION_HPP
