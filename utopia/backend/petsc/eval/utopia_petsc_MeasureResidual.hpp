#ifndef UTOPIA_PETSC_MEASURE_RESIDUAL_HPP
#define UTOPIA_PETSC_MEASURE_RESIDUAL_HPP

#include "utopia_petsc_Types.hpp"

#include "utopia_MeasureResidual.hpp"

namespace utopia {

    template <>
    class MeasureResidualComponents<PetscVector> final : public MeasureResidual<PetscVector> {
    public:
        using Super = utopia::MeasureResidual<PetscVector>;
        using Scalar = typename Traits<PetscVector>::Scalar;
        using SizeType = typename Traits<PetscVector>::SizeType;
        using ScalarArray = typename Traits<PetscVector>::ScalarArray;

        static std::string name() { return "MeasureResidualComponents"; }

        void read(Input &in) override;
        Scalar measure(const PetscVector &r) const override;

    private:
        int block_size_{1};
        int component_{0};
        bool verbose_{false};
    };

    template <>
    class MeasureResidualFactory<PetscVector> {
    public:
        static std::shared_ptr<MeasureResidual<PetscVector>> make(const std::string &type);
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_MEASURE_RESIDUAL_HPP
