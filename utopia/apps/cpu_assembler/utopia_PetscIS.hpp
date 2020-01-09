#ifndef UTOPIA_PETSC_IS_HPP
#define UTOPIA_PETSC_IS_HPP

#include "utopia_Traits.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    class PetscIS {
    public:
        using SizeType = Traits<PetscVector>::SizeType;
        using Destroy = std::function<PetscErrorCode(IS *)>;

        PetscIS(Destroy destroy_impl = ISDestroy)
        : is_(nullptr), destroy_impl_(destroy_impl), inds_(nullptr)
        {}

        PetscIS(IS is, Destroy destroy_impl = ISDestroy)
        : is_(is), destroy_impl_(destroy_impl), inds_(nullptr)
        {}

        const IS &raw_type() const
        {
            return is_;
        }

        IS &raw_type()
        {
            return is_;
        }

        inline void destroy()
        {
            if(is_) {
                destroy_impl_(&is_);
                is_ = nullptr;
            }
        }

        void read_lock()
        {
            if(inds_) return;
            ISGetIndices(is_, &inds_);
        }

        void read_unlock()
        {
            if(inds_) {
                ISRestoreIndices(is_, &inds_);
                inds_ = nullptr;
            }
        }

        inline SizeType size() const
        {
            SizeType s;
            ISGetSize(is_, &s);
            return s;
        }

        inline const SizeType &operator()(const SizeType &i) const {
            assert(inds_);
            return inds_[i];
        }

        inline const SizeType &operator[](const SizeType &i) const {
            assert(inds_);
            return inds_[i];
        }

    private:
        IS is_;
        Destroy destroy_impl_;
        const SizeType *inds_;
    };

    template<>
    class Traits<PetscIS> : public Traits<PetscVector> {};
}

#endif //UTOPIA_PETSC_IS_HPP
