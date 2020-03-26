#ifndef UTOPIA_FUNCTION_SPACE_VIEW_HPP
#define UTOPIA_FUNCTION_SPACE_VIEW_HPP

#include "utopia_Base.hpp"
#include "utopia_MemType.hpp"
#include "utopia_MeshView.hpp"
#include "utopia_Device.hpp"
#include "utopia_FunctionSpace.hpp"

namespace utopia {

    template<class MeshView_, int Dim>
    class FunctionSpaceView {
    public:
        using MeshView = MeshView_;

    private:
        MeshView mesh_;
    };

    template<class MeshView_>
    class FunctionSpaceView<MeshView_, 1> {
    public:
        using MeshView = MeshView_;
        using Elem     = typename MeshView::Elem;
        using SizeType = typename MeshView::SizeType;
        using DofIndex = utopia::ArrayView<SizeType, Elem::NNodes>;

        template<class Array>
        UTOPIA_INLINE_FUNCTION void dofs(const SizeType &element_idx, Array &indices) const
        {
            mesh_.nodes(element_idx, indices);
        }

        UTOPIA_INLINE_FUNCTION void elem(const SizeType &element_idx, Elem &elem) const
        {
            mesh_.elem(element_idx, elem);
        }

        FunctionSpaceView(const MeshView &mesh) : mesh_(mesh) {}

        UTOPIA_INLINE_FUNCTION const MeshView &mesh() const
        {
            return mesh_;
        }

        UTOPIA_INLINE_FUNCTION MeshView &mesh()
        {
            return mesh_;
        }

    private:
        MeshView mesh_;
    };

    template<class Elem_, class Comm, class ExecutionSpace_, typename...Args>
    class FunctionSpace< Mesh<Elem_, Comm, ExecutionSpace_, Uniform<Args...>>, 1> {
    public:
        using Elem           = Elem_;
        using Mesh           = utopia::Mesh<Elem, Comm, ExecutionSpace_, Uniform<Args...>>;
        using MeshView       = typename Mesh::ViewDevice;
        using ExecutionSpace = ExecutionSpace_;
        using SizeType       = typename Mesh::SizeType;
        using DofIndex       = utopia::ArrayView<std::size_t, Elem::NNodes>;
        using ViewDevice     = utopia::FunctionSpaceView<MeshView, 1>;
        using Device         = utopia::Device<TRILINOS>;
        using MemType        = utopia::Uniform<Args...>;

        template<class Fun>
        void each_element(Fun fun)
        {
            mesh_.each_element(fun);
        }

        ViewDevice view_device()
        {
            return ViewDevice(mesh_.view_device());
        }

        Range element_range() const
        {
            return mesh_.element_range();
        }

        FunctionSpace(const Mesh &mesh) : mesh_(mesh) {}
    private:
        Mesh mesh_;
    };

}

#endif //UTOPIA_FUNCTION_SPACE_VIEW_HPP
