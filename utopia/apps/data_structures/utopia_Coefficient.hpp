#ifndef UTOPIA_COEFFICIENT_HPP
#define UTOPIA_COEFFICIENT_HPP

#include "utopia_DeviceView.hpp"

// petsc
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    template <class FunctionSpaceView, class VectorView>
    class CoefficientView {
    public:
        using Elem = typename FunctionSpaceView::Elem;
        using SizeType = typename FunctionSpaceView::SizeType;

        template <class Values>
        void get(const Elem &e, Values &values) const {
            space_.local_coefficients(e, vec_, values);
        }

        template <class Values>
        void get(const Elem &e, const SizeType &var, Values &values) const {
            space_.local_coefficients(e, vec_, var, values);
        }

        CoefficientView(const FunctionSpaceView &space, const VectorView &vec) : space_(space), vec_(vec) {}

    private:
        FunctionSpaceView space_;
        VectorView vec_;
    };

    template <class FunctionSpace>
    class Coefficient
    // {};

    // template<class Elem, int Components>
    // class Coefficient<FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>>
    {
    public:
        using Mesh = typename FunctionSpace::Mesh;
        // using FunctionSpace = utopia::FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename FunctionSpace::Scalar;
        using ViewDevice =
            utopia::CoefficientView<typename FunctionSpace::ViewDevice, LocalViewDevice<const PetscVector, 1>>;

        Coefficient(const FunctionSpace &space, Vector &global_vector)
            : space_(space), local_vector_(std::make_shared<PetscVector>()) {
            space_.mesh().create_local_vector(*local_vector_);
            update(global_vector);
        }

        Coefficient(const FunctionSpace &space, const std::shared_ptr<PetscVector> &local_vector)
            : space_(space), local_vector_(local_vector) {}

        explicit Coefficient(const FunctionSpace &space) : space_(space), local_vector_(std::make_shared<Vector>()) {
            // FIXME
            space_.mesh().create_local_vector(*local_vector_);
        }

        void update(const PetscVector &global_vector) {
            // FIXME
            space_.mesh().global_to_local(global_vector, *local_vector_);
        }

        void set(const Scalar &val) { local_vector_->set(val); }

        ViewDevice view_device() const { return ViewDevice(space_.view_device(), *local_vector_); }

        inline const FunctionSpace &space() const { return space_; }

    private:
        const FunctionSpace &space_;
        std::shared_ptr<PetscVector> local_vector_;
    };

}  // namespace utopia

#endif  // UTOPIA_COEFFICIENT_HPP
