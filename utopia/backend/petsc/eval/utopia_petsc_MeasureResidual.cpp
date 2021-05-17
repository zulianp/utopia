#include "utopia_petsc_MeasureResidual.hpp"

#include "utopia_DeviceView.hpp"

namespace utopia {

    void MeasureResidualComponents<PetscVector>::read(Input &in) {
        in.get("block_size", block_size_);
        in.get("component", component_);
        in.get("verbose", verbose_);
    }

    Traits<PetscVector>::Scalar MeasureResidualComponents<PetscVector>::measure(const PetscVector &r) const {
        if (block_size_ == 1) {
            return Super::measure(r);
        }

        ScalarArray ret(block_size_);

        // FIXME not parallel and not on device
        {
            Read<PetscVector> read_r(r);
            auto rr = range(r);

            for (SizeType i = rr.begin(); i < rr.end(); i += block_size_) {
                for (int sub_i = 0; sub_i < block_size_; ++sub_i) {
                    auto x = r.get(i + sub_i);
                    auto x2 = x * x;
                    ret[sub_i] += x2;
                }
            }
        }

        r.comm().sum(block_size_, &ret[0]);

        for (int sub_i = 0; sub_i < block_size_; ++sub_i) {
            ret[sub_i] = std::sqrt(ret[sub_i]);
        }

        if (verbose_) {
            for (int i = 0; i < block_size_; ++i) {
                utopia::out() << "residual(" << i << "): " << ret[i] << '\n';
            }
        }

        return ret[component_];
    }

    std::shared_ptr<MeasureResidual<PetscVector>> MeasureResidualFactory<PetscVector>::make(const std::string &type) {
        if (type == MeasureResidualComponents<PetscVector>::name()) {
            return std::make_shared<MeasureResidualComponents<PetscVector>>();
        }

        return std::make_shared<MeasureResidual<PetscVector>>();
    }

}  // namespace utopia
