#ifndef UTOPIA_MARS_INTREPID2_DISCRETIZATION_HPP
#define UTOPIA_MARS_INTREPID2_DISCRETIZATION_HPP

#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_kokkos_Discretization.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_UniformFE.hpp"
#include "utopia_mars_FunctionSpace.hpp"

#include <memory>

namespace utopia {

    using mars_FS_t = utopia::mars::FunctionSpace;
    using mars_FE_t = utopia::kokkos::UniformFE<utopia::Traits<utopia::mars::FunctionSpace>::Scalar>;
    using mars_Discretization_t = Discretization<mars_FS_t, mars_FE_t>;

    template <>
    class Discretization<mars_FS_t, mars_FE_t> final : public utopia::kokkos::Discretization<mars_FS_t, mars_FE_t> {
    public:
        using FunctionSpace = mars_FS_t;
        using FE = mars_FE_t;
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
                                         const Part &part) override;

        void local_to_global(const std::vector<ScalarAccumulator> &acc,
                             std::vector<Scalar> &scalars,
                             const Part &part = all()) override;
        class Impl;

    private:
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_MARS_INTREPID2_DISCRETIZATION_HPP
