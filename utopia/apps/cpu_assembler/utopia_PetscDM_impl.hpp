#ifndef UTOPIA_PETSC_DM_IMPL_HPP
#define UTOPIA_PETSC_DM_IMPL_HPP


#include "utopia_PetscDM.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_Readable.hpp"
#include "utopia_PetscIS.hpp"
#include "utopia_Rename.hpp"
#include "utopia_petsc_Eval_Rename.hpp"
#include "utopia_petsc_Each.hpp"
#include "utopia_ArrayView.hpp"

#include <petscdm.h>
#include <petscdmda.h>

#include <csignal>

namespace utopia {

    template<int Dim>
    DM raw_type(const PetscDM<Dim> &dm);

    template<int Dim>
    class DMMirror  {
    public:
        static const std::size_t UDim = Dim;
        using SizeType = PetscInt;
        using IntArray = utopia::ArrayView<SizeType, UDim>;

        virtual ~DMMirror() {}

        static void dof_ownership_range(const DM &dm, SizeType &begin, SizeType &end)
        {
            Vec v;
            DMGetGlobalVector(dm, &v);
            VecGetOwnershipRange(v, &begin, &end);
            DMRestoreGlobalVector(dm, &v);
        }

        static SizeType get_dimension(DM dm)
        {
            SizeType ret;
            DMGetDimension(dm, &ret);
            return ret;
        }

        virtual void init(DM dm)
        {
            dof_ownership_range(dm, dof_range_begin, dof_range_end);
        }

        SizeType dof_range_begin;
        SizeType dof_range_end;
    };

    template<int Dim>
    class DMDAMirror final : public DMMirror<Dim> {
    public:
        using Super = utopia::DMMirror<Dim>;
        static const std::size_t UDim = Dim;
        using SizeType = PetscInt;
        using IntArray = utopia::ArrayView<SizeType, UDim>;
        using Point    = typename PetscElem<Dim>::Point;

        static void get_dims(DM dm, IntArray &dims)
        {
            utopia::ArrayView<SizeType, 3> dims_buff;
            auto ierr = DMDAGetInfo(dm, nullptr, &dims_buff[0], &dims_buff[1], &dims_buff[2], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);

            for(int d = 0; d < Dim; ++d) {
                dims[d] = dims_buff[d];
            }
        }

        static void get_corners(DM dm, IntArray &start, IntArray &extent)
        {
            utopia::ArrayView<SizeType, 3> start_buff;
            utopia::ArrayView<SizeType, 3> extent_buff;

            DMDAGetCorners(dm,
                &start_buff[0],  &start_buff[1],  &start_buff[2],
                &extent_buff[0], &extent_buff[1], &extent_buff[2]
            );

            for(int d = 0; d < Dim; ++d) {
                start[d]  = start_buff[d];
                extent[d] = extent_buff[d];
            }
        }

        static void get_ghost_corners(DM dm, IntArray &start, IntArray &extent)
        {
            utopia::ArrayView<SizeType, 3> start_buff;
            utopia::ArrayView<SizeType, 3> extent_buff;

            DMDAGetGhostCorners(dm,
                &start_buff[0],  &start_buff[1],  &start_buff[2],
                &extent_buff[0], &extent_buff[1], &extent_buff[2]
            );

            for(int d = 0; d < Dim; ++d) {
                start[d]  = start_buff[d];
                extent[d] = extent_buff[d];
            }
        }

        static SizeType get_dof(DM dm)
        {
            SizeType ret;
            DMDAGetDof(dm, &ret);
            return ret;
        }

        void init(DM dm) override
        {
            Super::init(dm);
            get_dims(dm, dims);
            get_corners(dm, corners_begin, corners_extent);
            get_ghost_corners(dm, ghost_corners_begin, ghost_corners_extent);
            n_components = get_dof(dm);

            describe();
        }

        DMDAMirror()
        {
            box_min.set(0.0);
            box_max.set(1.0);
        }

        // void describe(std::ostream &os = std::cout)
        void describe() const
        {
            disp("dims");
            disp(dims);

            disp("corners_begin");
            disp(corners_begin);

            disp("corners_extent");
            disp(corners_extent);


            disp("ghost_corners_begin");
            disp(ghost_corners_begin);

            disp("ghost_corners_extent");
            disp(ghost_corners_extent);

            disp("n_components");
            disp(n_components);

            disp("dof_range_begin");
            disp(this->dof_range_begin);

            disp("dof_range_end");
            disp(this->dof_range_end);

            //geometry
            disp("box_min");
            disp(box_min);

            disp("box_max");
            disp(box_max);
        }

        IntArray dims;
        IntArray corners_begin;
        IntArray corners_extent;

        IntArray ghost_corners_begin;
        IntArray ghost_corners_extent;

        SizeType n_components;


        //geometry
        Point box_min, box_max;
    };

    template<int Dim>
    class DMDAElements {
    public:
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        using Scalar   = utopia::Traits<PetscVector>::Scalar;

        DMDAElements(const PetscDM<Dim> &dm_impl)
        : dm(raw_type(dm_impl))
        {
            DMDAGetElements(dm,&ne,&nc,&e);
            e_global.resize(ne*nc);

            ISLocalToGlobalMapping map;
            DMGetLocalToGlobalMapping(dm,&map);
            ISLocalToGlobalMappingApplyBlock(map, ne*nc, e, &e_global[0]);
        }

        ~DMDAElements()
        {
            DMDARestoreElements(dm,&ne,&nc,&e);
        }

        DM dm;
        SizeType ne,nc;
        const SizeType   *e;
        std::vector<SizeType> e_global;
    };

    template<int Dim>
    class DMDANodes {
    public:
        using SizeType = PetscInt;
        using Scalar   = PetscScalar;
        using PetscNode = utopia::PetscNode<Dim>;
        using Point = typename utopia::PetscDM<Dim>::Point;

        template<class Array>
        void dims(Array &arr) const
        {
            auto ierr = DMDAGetInfo(dm, nullptr, &arr[0], &arr[1], &arr[2], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);
        }

        template<class Fun>
        void each_with_ghosts(Fun f)
        {
            SizeType dims[3], start[3], extent[3];
            DMDAGetGhostCorners(dm,
                &start[0],  &start[1],  &start[2],
                &extent[0], &extent[1], &extent[2]
                );



            this->dims(dims);

            each(dims, start, extent, f);
        }

        template<class Fun>
        void each(
            const SizeType dims[3],
            const SizeType start[3],
            const SizeType extent[3], Fun f) const
        {
            switch(Dim) {
                case 1:
                {
                    for(SizeType i = 0; i < extent[0]; ++i) {
                        SizeType idx = i + start[0];
                        PetscNode node(*this, idx);
                        f(node);
                    }

                    break;
                }

                case 2:
                {
                    for(SizeType j = 0; j < extent[1]; ++j) {
                        for(SizeType i = 0; i < extent[0]; ++i) {
                            SizeType idx = (j + start[1]) * dims[0] + i + start[0];
                            PetscNode node(*this, idx);
                            f(node);
                        }
                    }
                    break;
                }

                case 3:
                {
                    for(SizeType k = 0; k < extent[2]; ++k) {
                        const SizeType k_offset = (k + start[2]) * dims[2];

                        for(SizeType j = 0; j < extent[1]; ++j) {
                            const SizeType j_offset = (j + start[1]) * dims[0];

                            for(SizeType i = 0; i < extent[0]; ++i) {
                                SizeType idx = k_offset + j_offset + i + start[0];

                                PetscNode node(*this, idx);
                                f(node);
                            }
                        }
                    }
                    break;
                }
                default:
                {
                    assert(false);
                    break;
                }
            }
        }

        inline bool is_ghost(const SizeType global_idx) const
        {
            return !range.inside(global_idx);
        }

        inline SizeType dim() const
        {
            SizeType ret;
            DMGetDimension(dm, &ret);
            return ret;
        }

        DMDANodes(const PetscDM<Dim> &wrapper)
        : dm(raw_type(wrapper)),
          range(wrapper.mirror().dof_range_begin, wrapper.mirror().dof_range_end)
        {}

    private:
        DM dm;
        Range range;
    };

    template<int Dim>
    class PetscDM<Dim>::Impl {
    public:
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        using Scalar   = utopia::Traits<PetscVector>::Scalar;
        using Point    = typename PetscElem<Dim>::Point;

        Impl(const PetscCommunicator &comm)
        : comm(comm), dm(nullptr)
        {}

        ~Impl() {
            destroy();
        }

        void init(
            const PetscCommunicator &comm,
            const std::array<SizeType, UDim> &arr,
            DMDAStencilType stencil_type,
            const SizeType n_components)
        {
            destroy();
            this->comm = comm;

            const SizeType dim = arr.size();

            switch(dim)
            {
                case 1:
                {
                    DMDACreate1d(comm.get(),
                        DM_BOUNDARY_NONE,
                        arr[0],
                            n_components, //dofs per node
                            1, //stencil width
                            nullptr,
                            &dm
                            );

                    break;
                }

                case 2:
                {
                    DMDACreate2d(comm.get(),
                        DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
                        stencil_type,
                        arr[0], arr[1],
                        PETSC_DECIDE, PETSC_DECIDE,
                            n_components, //dofs per node
                            1, //stencil width
                            nullptr,
                            nullptr,
                            &dm
                            );

                    break;
                }

                case 3:
                {
                    DMDACreate3d(comm.get(),
                        DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
                        stencil_type,
                        arr[0], arr[1], arr[2],
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                            n_components, //dofs per node
                            1, //stencil width
                            nullptr,
                            nullptr,
                            nullptr,
                            &dm
                            );

                    break;
                }

                default:
                {
                    assert(false);
                }
            }
        }

        void init_uniform(
            const PetscCommunicator &comm,
            const std::array<SizeType, UDim> &arr,
            const std::array<Scalar, UDim> &box_min,
            const std::array<Scalar, UDim> &box_max,
            const SizeType &n_components,
            DMDAElementType elem_type    = DMDA_ELEMENT_Q1,
            DMDAStencilType stencil_type = DMDA_STENCIL_BOX)
        {
            init(comm, arr, stencil_type, n_components);

            Scalar min_x = 0, min_y = 0, min_z = 0;
            Scalar max_x = 1, max_y = 1, max_z = 1;

            min_x = box_min[0];
            max_x = box_max[0];

            for(int d = 0; d < Dim; ++d) {
                this->mirror.box_min[d] = box_min[d];
                this->mirror.box_max[d] = box_max[d];
            }

            const auto d = arr.size();
            if(d > 1) {
                min_y = box_min[1];
                max_y = box_max[1];
            }

            if(d > 2) {
                min_z = box_min[2];
                max_z = box_max[2];
            }


            DMDASetElementType(dm, elem_type);
            DMSetUp(dm);
            DMDASetUniformCoordinates(dm, min_x, max_x, min_y, max_y, min_z, max_z);
        }

        void destroy()
        {
            if(dm) {
                elements = nullptr;
                nodes = nullptr;
                DMDestroy(&dm);
                dm = nullptr;
            }
        }

        template<class Array>
        void local_node_grid_coord_no_ghost(const SizeType &idx, Array &tensor_index) const
        {
            SizeType current = idx;
            for(int i = 0; i < Dim; ++i) {
                const int next = current / mirror.corners_extent[i];
                tensor_index[i] = current - next * mirror.corners_extent[i] + mirror.corners_begin[i];
                current = next;
            }
        }

        template<class Array>
        void local_node_grid_coord(const SizeType &idx, Array &tensor_index) const
        {
            SizeType current = idx;
            for(int i = 0; i < Dim; ++i) {
                const int next = current / mirror.ghost_corners_extent[i];
                tensor_index[i] = current - next * mirror.ghost_corners_extent[i];
                current = next;
            }
        }

        void local_node_idx_coord(const SizeType &idx, Point &p)
        {
            SizeType tensor_index[Dim];

            SizeType current = idx;
            for(int i = 0; i < Dim; ++i) {
                const int next = current / mirror.ghost_corners_extent[i];
                tensor_index[i] = current - next * mirror.ghost_corners_extent[i];
                current = next;
            }

            for(int i = 0; i < Dim; ++i) {
                SizeType c = tensor_index[i] + mirror.ghost_corners_begin[i];
                p[i] = c * (mirror.box_max[i] - mirror.box_min[i])/(mirror.dims[i] - 1) + mirror.box_min[i];
            }
        }

        PetscCommunicator comm;
        DM dm;

        std::unique_ptr<DMDAElements<Dim>> elements;
        std::unique_ptr<DMDANodes<Dim>> nodes;
        DMDAMirror<Dim> mirror;
    };

    template<int Dim>
    PetscDM<Dim>::PetscDM()
    : impl_(utopia::make_unique<Impl>(PetscCommunicator()))
    {}

    template<int Dim>
    DM raw_type(const PetscDM<Dim> &dm)
    {
        return dm.impl().dm;
    }

    template<int Dim>
    Range local_node_range(const PetscDM<Dim> &dm)
    {
        return dm.local_node_range();
    }

    template<int Dim>
    PetscDM<Dim>::PetscDM(
        const PetscCommunicator &comm,
        const std::array<SizeType, UDim> &dims,
        const std::array<Scalar, UDim> &box_min,
        const std::array<Scalar, UDim> &box_max,
        const SizeType &n_components)
    {
        build(comm, dims, box_min, box_max, n_components);
    }

    template<int Dim>
    void PetscDM<Dim>::build(
        const PetscCommunicator &comm,
        const std::array<SizeType, UDim> &dims,
        const std::array<Scalar, UDim> &box_min,
        const std::array<Scalar, UDim> &box_max,
        const SizeType &n_components)
    {
        impl_ = utopia::make_unique<Impl>(comm);
        impl_->init_uniform(comm, dims, box_min, box_max, n_components);

        //FIXME move into init uniform
        impl_->mirror.init(impl_->dm);

        impl_->elements = utopia::make_unique<DMDAElements<Dim>>(*this);
        impl_->nodes    = utopia::make_unique<DMDANodes<Dim>>(*this);
    }

    template<int Dim>
    PetscDM<Dim>::~PetscDM()
    {}

    template<int Dim>
    Range PetscDM<Dim>::local_element_range() const
    {
        return Range(0, impl_->elements->ne);
    }

    template<int Dim>
    void PetscDM<Dim>::elem(const SizeType &idx, Elem &e) const
    {
        this->nodes(idx, e.nodes());
        e.idx(idx);
    }

    template<int Dim>
    void PetscDM<Dim>::nodes(const SizeType &idx, NodeIndex &nodes) const
    {
        assert(idx < impl_->elements->ne);
        assert(idx >= 0);
        assert(!impl_->elements->e_global.empty());
        nodes = NodeIndex(&impl_->elements->e_global[idx*impl_->elements->nc], impl_->elements->nc);
    }

    template<int Dim>
    void PetscDM<Dim>::nodes_local(const SizeType &idx, NodeIndex &nodes) const
    {
        assert(idx < impl_->elements->ne);
        assert(idx >= 0);
        assert(impl_->elements->e);
        nodes = NodeIndex(&impl_->elements->e[idx*impl_->elements->nc], impl_->elements->nc);
    }

    template<int Dim>
    Range PetscDM<Dim>::local_node_range() const
    {
        return Range(mirror().dof_range_begin, mirror().dof_range_end);
    }

    template<int Dim>
    void PetscDM<Dim>::create_matrix(PetscMatrix &mat) const
    {
        mat.destroy();
        DMCreateMatrix(impl_->dm, &mat.raw_type());
    }

    template<int Dim>
    void PetscDM<Dim>::create_vector(PetscVector &vec) const
    {
        vec.destroy();
        DMCreateGlobalVector(impl_->dm, &vec.raw_type());
    }

    template<int Dim>
    void PetscDM<Dim>::local_to_global(const PetscVector &local,  PetscVector &global) const
    {
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0) //DMA-INCOMPLETE
        DMLocalToGlobal(impl_->dm, local.raw_type(), ADD_VALUES, global.raw_type());
#else
        DMLocalToGlobalBegin(impl_->dm, local.raw_type(), ADD_VALUES, global.raw_type());
        DMLocalToGlobalEnd(impl_->dm, local.raw_type(), ADD_VALUES, global.raw_type());
#endif
    }

    template<int Dim>
    void PetscDM<Dim>::global_to_local(const PetscVector &global, PetscVector &local) const
    {
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 0) //DMA-INCOMPLETE
        DMGlobalToLocal(impl_->dm, global.raw_type(), INSERT_VALUES, local.raw_type());
#else
        DMGlobalToLocalBegin(impl_->dm, global.raw_type(), INSERT_VALUES, local.raw_type());
        DMGlobalToLocalEnd(impl_->dm, global.raw_type(), INSERT_VALUES, local.raw_type());
#endif
    }

    template<int Dim>
    void PetscDM<Dim>::describe() const
    {
        PetscErrorCode err = 0;
        SizeType dim, M, N, P, m, n, p, dof, s;
        DMDAStencilType st;
        err = DMDAGetInfo(impl_->dm, &dim, &M, &N, &P, &m, &n, &p, &dof, &s, nullptr, nullptr, nullptr, &st); assert(err == 0);

        int size = impl_->comm.size();
        int rank = impl_->comm.rank();

        impl_->comm.barrier();

        Vec v;
        DMGetGlobalVector(impl_->dm, &v);

        PetscInt n_b, n_e;
        VecGetOwnershipRange(v, &n_b, &n_e);
        Range range(n_b, n_e);

        for(int i = 0; i < size; ++i) {
            if(i == rank) {
                std::cout << "--------------------------------------------\n";
                std::cout << i << ")\ndim: " << dim << "\nM: " << M << "\nN: " << N << "\nP: " << P << "\nm: " << m << "\nn: " << n << "\np: " << p << "\ndof: " << dof << "\ns: " << s <<  std::endl;

                DMDANodes<Dim> nodes(*this);
                nodes.each_with_ghosts([range](const Node &node) {

                    std::cout << node.idx() << " ";
                    if(range.inside(node.idx())) {

                    } else {
                        std::cout << "g";
                    }

                    std::cout << std::endl;
                });

                std::cout << "--------------------------------------------\n";
            }

            impl_->comm.barrier();
        }

        DMRestoreGlobalVector(impl_->dm, &v);
    }

    static void print_corner(
        const std::vector<SizeType> &c,
        std::ostream &os = std::cout)
    {
        const SizeType dim = c.size();
        for(SizeType d = 0; d < dim; ++d) {
            os << c[d] << " ";
        }

        os << std::endl;
    }

    template<int Dim>
    void PetscDM<Dim>::create_local_vector(PetscVector &vec) const
    {
        vec.destroy();
        auto err = DMCreateLocalVector(impl_->dm, &raw_type(vec)); assert(err == 0);
    }

    template<int Dim>
    bool PetscNode<Dim>::is_ghost() const
    {
        return nodes_.is_ghost(idx());
    }

    template<int Dim>
    bool PetscDM<Dim>::is_ghost(const SizeType &global_node_idx) const
    {
        return impl_->nodes->is_ghost(global_node_idx);
    }

    template<int Dim>
    bool PetscDM<Dim>::is_local_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const
    {
        std::array<SizeType, 3> tensor_index = {0, 0, 0};
        impl_->local_node_grid_coord_no_ghost(idx, tensor_index);

        const auto &mirror = impl_->mirror;

        switch(b_id) {
            case SideSet::left():
            {
                return tensor_index[0] == 0;
            }

            case SideSet::right():
            {
                return tensor_index[0] == (mirror.dims[0] - 1);
            }

            case SideSet::bottom():
            {
                return tensor_index[1] == 0;
            }

            case SideSet::top():
            {
                return tensor_index[1] == (mirror.dims[1] - 1);
            }

            case SideSet::back():
            {
                return tensor_index[2] == 0;
            }

            case SideSet::front():
            {
                return tensor_index[2] == (mirror.dims[2] - 1);
            }

            default:
            {
                break;
            }
        }

        return false;
    }

    template<int Dim>
    void PetscDM<Dim>::cell_point(const SizeType &idx, Point &translation)
    {
        SizeType v0 = impl_->elements->e[idx*impl_->elements->nc];
        impl_->local_node_idx_coord(v0, translation);
    }

    template<int Dim>
    void PetscDM<Dim>::cell_size(const SizeType &, Point &cell_size)
    {
        const auto &mirror = impl_->mirror;
        for(int d = 0; d < Dim; ++d) {
            cell_size[d] = (impl_->mirror.box_max[d] - impl_->mirror.box_min[d])/(mirror.dims[d] - 1);
        }
    }

    template<int Dim>
    PetscCommunicator &PetscDM<Dim>::comm()
    {
        return impl_->comm;
    }

    template<int Dim>
    const PetscCommunicator &PetscDM<Dim>::comm() const
    {
        return impl_->comm;
    }

    template<int Dim>
    typename PetscDM<Dim>::Scalar PetscDM<Dim>::min_spacing() const
    {
        const auto &mirror = impl_->mirror;
        Scalar min_h = (impl_->mirror.box_max[0] - mirror.box_min[0])/mirror.dims[0];
        for(int i = 1; i < Dim; ++i) {
            Scalar h = (mirror.box_max[i] - mirror.box_min[i])/mirror.dims[i];
            min_h = std::min(min_h, h);
        }

        return min_h;
    }

    template<int Dim>
    const DMDAMirror<Dim> &PetscDM<Dim>::mirror() const
    {
        return impl_->mirror;
    }

    template<int Dim>
    typename PetscDM<Dim>::SizeType PetscDM<Dim>::n_nodes() const
    {
        const auto &mirror = impl_->mirror;

        SizeType ret = mirror.dims[0];
        for(int i = 1; i < Dim; ++i) {
            ret *= mirror.dims[i];
        }

        return ret;
    }

    template<int Dim>
    typename PetscDM<Dim>::SizeType PetscDM<Dim>::n_local_nodes() const
    {
        return impl_->mirror.dof_range_end - impl_->mirror.dof_range_begin;
    }

    template<int Dim>
    typename PetscDM<Dim>::SizeType PetscDM<Dim>::n_elements() const
    {
        const auto &mirror = impl_->mirror;

        SizeType ret = mirror.dims[0] - 1;

        for(int i = 1; i < Dim; ++i) {
            ret *= (mirror.dims[i] - 1);
        }

        return ret;
    }

    template<int Dim>
    typename PetscDM<Dim>::SizeType PetscDM<Dim>::n_components() const
    {
        return impl_->mirror.n_components;
    }

    template<int Dim>
    void PetscDM<Dim>::set_field_name(const SizeType &nf, const std::string &name)
    {
        DMDASetFieldName(impl_->dm, nf, name.c_str());
    }

}

#endif //UTOPIA_PETSC_DM_IMPL_HPP
