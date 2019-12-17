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

#include <petscdm.h>
#include <petscdmda.h>

#include <csignal>

namespace utopia {


    //DMDAGetGlobalIndices for local 2 global mapping including ghost nodes
    //https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex42.c.html
    //https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/FE/index.html

    template<int Dim>
    DM raw_type(const PetscDM<Dim> &dm);

    template<int Dim>
    class PetscDMImpl {
    public:
        using SizeType = Traits<PetscVector>::SizeType;
        using Scalar   = Traits<PetscVector>::Scalar;

        template<class Array>
        inline static void dims(const DM &dm, Array &arr)
        {
            auto ierr = DMDAGetInfo(dm, nullptr, &arr[0], &arr[1], &arr[2], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);
        }

        static Range local_node_range(const DM &dm)
        {
            Vec v;
            DMGetGlobalVector(dm, &v);

            PetscInt n_b, n_e;
            VecGetOwnershipRange(v, &n_b, &n_e);
            Range range;
            range.set(n_b, n_e);

            DMRestoreGlobalVector(dm, &v);
            return range;
        }

        static void local_node_ranges(const DM &dm, SizeType *begin, SizeType *end)
        {
            DMDAGetCorners(dm,
                &begin[0],  &begin[1],  &begin[2],
                &end[0],    &end[1],    &end[2]
                );

            for(int d = 0; d < 3; ++d) {
                end[d] += begin[d];
            }
        }

        static void local_element_ranges(const DM &dm, SizeType *begin, SizeType *end)
        {
            DMDAGetElementsCorners(dm, &begin[0],  &begin[1],  &begin[2]);
            DMDAGetElementsSizes(dm,   &end[0],    &end[1],    &end[2]);

            for(int d = 0; d < 3; ++d) {
                end[d] += begin[d];
            }
        }

        // static void box(const DM &dm, Scalar *min, Scalar *max)
        // {
        //     DMDAGetBoundingBox(dm, min, max);
        // }

        static SizeType n_local_nodes_with_ghosts(const DM &dm)
        {
            SizeType start[3], extent[3];
            DMDAGetGhostCorners(dm,
                &start[0],  &start[1],  &start[2],
                &extent[0], &extent[1], &extent[2]
                );

            SizeType ret = extent[0];
            SizeType d = dim(dm);

            for(SizeType i = 1; i < d; ++i) {
                ret *= extent[i];
            }

            return ret;
        }

        inline static SizeType dim(const DM &dm)
        {
            SizeType ret;
            DMGetDimension(dm, &ret);
            return ret;
        }

        template<class Array>
        inline static void point(const DM &dm, const SizeType &i, const SizeType &j, const SizeType &k, Array &array)
        {
            assert( dim(dm) == 3 );
            DMDAGetCellPoint(dm, i, j, k, &array[0]);
        }

        template<class Array>
        inline static void point(const DM &dm, const SizeType &i, const SizeType &j, Array &array)
        {
            assert( dim(dm) == 2 );
            DMDAGetCellPoint(dm, i, j, 0, &array[0]);
        }

        template<class Array>
        inline static void point(const DM &dm, const SizeType &i, Array &array)
        {
            assert( dim(dm) == 1 );
            DMDAGetCellPoint(dm, i, 0, 0, &array[0]);
        }
    };

    template<int Dim>
    class PetscDMElements {
    public:
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        using Scalar   = utopia::Traits<PetscVector>::Scalar;

        PetscDMElements(const PetscDM<Dim> &dm_impl)
        : dm(raw_type(dm_impl))
        {
            DMDAGetElements(dm,&ne,&nc,&e);
            e_global.resize(ne*nc);

            ISLocalToGlobalMapping map;
            DMGetLocalToGlobalMapping(dm,&map);
            ISLocalToGlobalMappingApplyBlock(map, ne*nc, e, &e_global[0]);
        }

        template<class Fun>
        void each(Fun fun) const
        {
            for(typename PetscDM<Dim>::SizeType i = 0; i < ne; i++) {
                // PetscDM<Dim>::Elem e(*this, i);
                // fun(e);
            }
        }

        ~PetscDMElements()
        {
            DMDARestoreElements(dm,&ne,&nc,&e);
        }

        DM dm;
        SizeType ne,nc;
        const SizeType   *e;
        std::vector<SizeType> e_global;
    };

    template<int Dim>
    class PetscDMNodes {
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
        void each(Fun f)
        {
            SizeType dims[3], start[3], extent[3];

            this->dims(dims);

            DMDAGetCorners(dm,
                &start[0],  &start[1],  &start[2],
                &extent[0], &extent[1], &extent[2]
                );

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

        PetscDMNodes(const PetscDM<Dim> &impl)
        : dm(raw_type(impl)), range(PetscDMImpl<Dim>::local_node_range(dm))
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
            DMDAStencilType stencil_type)
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
                            1, //dofs per node
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
                            1, //dofs per node
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
                            1, //dofs per node
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
            DMDAElementType elem_type    = DMDA_ELEMENT_Q1,
            DMDAStencilType stencil_type = DMDA_STENCIL_BOX)
        {
            init(comm, arr, stencil_type);

            Scalar min_x = 0, min_y = 0, min_z = 0;
            Scalar max_x = 1, max_y = 1, max_z = 1;

            min_x = box_min[0];
            max_x = box_max[0];

            for(int d = 0; d < Dim; ++d) {
                this->box_min_[d] = box_min[d];
                this->box_max_[d] = box_max[d];
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
        void global_node_grid_coord(const SizeType &idx, Array &tensor_index) const
        {
            SizeType dims[3];
            PetscDMImpl<Dim>::dims(dm, dims);

            SizeType current = idx;
            for(int i = 0; i < Dim; ++i) {
                const int next = current / dims[i];
                tensor_index[i] = current - next * dims[i];
                current = next;
            }
        }

        template<class Array>
        void local_node_grid_coord_no_ghost(const SizeType &idx, Array &tensor_index) const
        {
            SizeType dims[3], start[3], extent[3];
            DMDAGetCorners(dm,
                &start[0],  &start[1],  &start[2],
                &extent[0], &extent[1], &extent[2]
            );

            PetscDMImpl<Dim>::dims(dm, dims);

            SizeType current = idx;
            for(int i = 0; i < Dim; ++i) {
                const int next = current / extent[i];
                tensor_index[i] = current - next * extent[i] + start[i];
                current = next;
            }
        }

        template<class Array>
        void local_node_grid_coord(const SizeType &idx, Array &tensor_index) const
        {
            SizeType dims[3], start[3], extent[3];
            DMDAGetGhostCorners(dm,
                &start[0],  &start[1],  &start[2],
                &extent[0], &extent[1], &extent[2]
            );

            PetscDMImpl<Dim>::dims(dm, dims);

            SizeType current = idx;
            for(int i = 0; i < Dim; ++i) {
                const int next = current / extent[i];
                tensor_index[i] = current - next * extent[i];
                current = next;
            }
        }

        void local_node_idx_coord(const SizeType &idx, Point &p)
        {
            SizeType dims[3], start[3], extent[3], tensor_index[Dim];
            DMDAGetGhostCorners(dm,
                &start[0],  &start[1],  &start[2],
                &extent[0], &extent[1], &extent[2]
            );

            PetscDMImpl<Dim>::dims(dm, dims);

            SizeType current = idx;
            for(int i = 0; i < Dim; ++i) {
                const int next = current / extent[i];
                tensor_index[i] = current - next * extent[i];
                current = next;
            }

            for(int i = 0; i < Dim; ++i) {
                SizeType c = tensor_index[i] + start[i];
                p[i] = c * (box_max_[i] - box_min_[i])/(dims[i] - 1) + box_min_[i];
            }
        }

        PetscCommunicator comm;
        DM dm;

        std::unique_ptr<PetscDMElements<Dim>> elements;
        std::unique_ptr<PetscDMNodes<Dim>> nodes;
        Point box_min_, box_max_;
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
        const std::array<Scalar, UDim> &box_max)
    : impl_(utopia::make_unique<Impl>(comm))
    {
        impl_->init_uniform(comm, dims, box_min, box_max);
        impl_->elements = utopia::make_unique<PetscDMElements<Dim>>(*this);
        impl_->nodes = utopia::make_unique<PetscDMNodes<Dim>>(*this);
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
    }

    template<int Dim>
    void PetscDM<Dim>::nodes(const SizeType &idx, NodeIndex &nodes) const
    {
        nodes = NodeIndex(&impl_->elements->e_global[idx*impl_->elements->nc], impl_->elements->nc);
    }


    template<int Dim>
    void PetscDM<Dim>::nodes_local(const SizeType &idx, NodeIndex &nodes) const
    {
        nodes = NodeIndex(&impl_->elements->e[idx*impl_->elements->nc], impl_->elements->nc);
    }

    template<int Dim>
    void PetscDM<Dim>::each_element(const std::function<void(const Elem &)> &f)
    {
        impl_->elements->each(f);
    }

    template<int Dim>
    void PetscDM<Dim>::each_node(const std::function<void(const Node &)> &f)
    {
        impl_->nodes->each(f);
    }

    template<int Dim>
    void PetscDM<Dim>::each_node_with_ghosts(const std::function<void (const Node &)> &f)
    {
        impl_->nodes->each_with_ghosts(f);
    }

    template<int Dim>
    Range PetscDM<Dim>::local_node_range() const
    {
        return PetscDMImpl<Dim>::local_node_range(impl_->dm);
    }

    template<int Dim>
    void PetscDM<Dim>::box(Scalar *min, Scalar *max) const
    {
        for(int i = 0; i < Dim; ++i) {
            min[i] = impl_->box_min_[i];
            max[i] = impl_->box_max_[i];
        }
    }

    template<int Dim>
    void PetscDM<Dim>::local_node_ranges(SizeType *begin, SizeType *end) const
    {
        PetscDMImpl<Dim>::local_node_ranges(impl_->dm, begin, end);
    }

    template<int Dim>
    void PetscDM<Dim>::local_element_ranges(SizeType *begin, SizeType *end) const
    {
        PetscDMImpl<Dim>::local_element_ranges(impl_->dm, begin, end);
    }

    template<int Dim>
    typename PetscDM<Dim>::SizeType PetscDM<Dim>::n_local_nodes_with_ghosts() const
    {
        return PetscDMImpl<Dim>::n_local_nodes_with_ghosts(impl_->dm);
    }

    template<int Dim>
    void PetscDM<Dim>::create_matrix(PetscMatrix &mat)
    {
        mat.destroy();
        DMCreateMatrix(impl_->dm, &mat.raw_type());
    }

    template<int Dim>
    void PetscDM<Dim>::create_vector(PetscVector &vec)
    {
        vec.destroy();
        DMCreateGlobalVector(impl_->dm, &vec.raw_type());
    }

    template<int Dim>
    void PetscDM<Dim>::dims(std::array<SizeType, UDim> &arr) const
    {
        switch(Dim) {
            case 1:
            {
                auto ierr = DMDAGetInfo(impl_->dm, nullptr, &arr[0], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);
                break;
            }

            case 2:
            {
                auto ierr = DMDAGetInfo(impl_->dm, nullptr, &arr[0], &arr[1], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);
                break;
            }

            case 3:
            {
                auto ierr = DMDAGetInfo(impl_->dm, nullptr, &arr[0], &arr[1], &arr[2], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);
                break;
            }

            default: {
                static_assert(UDim >= 1 && UDim <= 3, "valid mesh sizes are 1, 2, 3");
                break;
            }
        }

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

                PetscDMNodes<Dim> nodes(*this);
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

    // template<int Dim>
    // void PetscDM<Dim>::create_local_matrix(PetscMatrix &mat)
    // {

    // }

    template<int Dim>
    void PetscDM<Dim>::create_local_vector(PetscVector &vec)
    {
        vec.destroy();
        auto err = DMCreateLocalVector(impl_->dm, &raw_type(vec)); assert(err == 0);
    }

    // template<int Dim>
    // void PetscDM<Dim>::local_to_global(PetscMatrix &local, PetscMatrix &global)
    // {

    // }

    template<int Dim>
    void PetscDM<Dim>::local_to_global(PetscVector &local, PetscVector &global)
    {
        auto err = DMLocalToGlobal(impl_->dm, raw_type(local), ADD_VALUES, raw_type(global)); assert(err == 0);
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
    bool PetscDM<Dim>::is_boundary(const SizeType &global_node_idx) const
    {
        std::array<SizeType, Dim> tensor_index;
        SizeType dims[3];
        impl_->global_node_grid_coord(global_node_idx, tensor_index);
        PetscDMImpl<Dim>::dims(impl_->dm, dims);

        for(int d = 0; d < Dim; ++d) {
            if(tensor_index[d] == 0 || tensor_index[d] == dims[d] -1) {
                return true;
            }
        }

        return false;
    }

    template<int Dim>
    bool PetscDM<Dim>::is_local_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const
    {
        std::array<SizeType, Dim> tensor_index;
        SizeType dims[3];
        impl_->local_node_grid_coord_no_ghost(idx, tensor_index);
        PetscDMImpl<Dim>::dims(impl_->dm, dims);

        switch(b_id) {
            case SideSet::left():
            {
                return tensor_index[0] == 0;
            }

            case SideSet::right():
            {
                return tensor_index[0] == (dims[0] - 1);
            }

            case SideSet::bottom():
            {
                return tensor_index[1] == 0;
            }

            case SideSet::top():
            {
                return tensor_index[1] == (dims[1] - 1);
            }

            case SideSet::back():
            {
                return tensor_index[2] == 0;
            }

            case SideSet::front():
            {
                return tensor_index[2] == (dims[2] - 1);
            }

            default:
            {
                break;
            }
        }

        return false;

    }

    template<int Dim>
    bool PetscDM<Dim>::is_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const
    {
        assert(false && "FIXME");

        std::array<SizeType, Dim> tensor_index;
        SizeType dims[3];
        impl_->global_node_grid_coord(idx, tensor_index);
        PetscDMImpl<Dim>::dims(impl_->dm, dims);

        switch(b_id) {
            case SideSet::left():
            {
                return tensor_index[0] == 0;
            }

            case SideSet::right():
            {
                return tensor_index[0] == (dims[0] - 1);
            }

            case SideSet::bottom():
            {
                return tensor_index[1] == 0;
            }

            case SideSet::top():
            {
                return tensor_index[1] == (dims[1] - 1);
            }

            case SideSet::back():
            {
                return tensor_index[2] == 0;
            }

            case SideSet::front():
            {
                return tensor_index[2] == (dims[2] - 1);
            }

            default:
            {
                break;
            }
        }

        return false;
    }

    template<int Dim>
    SideSet::BoundaryIdType PetscDM<Dim>::boundary_id(const SizeType &global_node_idx) const
    {
        std::array<SizeType, Dim> tensor_index;
        SizeType dims[3];
        impl_->global_node_grid_coord(global_node_idx, tensor_index);
        PetscDMImpl<Dim>::dims(impl_->dm, dims);

        if(tensor_index[0] == 0) {
            return SideSet::left();
        }

        if(tensor_index[0] == (dims[0] - 1)) {
            return SideSet::right();
        }

        if(Dim > 1) {
            if(tensor_index[1] == 0) {
                return SideSet::bottom();
            }

            if(tensor_index[1] == (dims[1] - 1)) {
                return SideSet::top();
            }

            if(Dim > 2) {
                if(tensor_index[2] == 0) {
                    return SideSet::back();
                }

                if(tensor_index[2] == (dims[2] - 1)) {
                    return SideSet::front();
                }
            }
        }

        return SideSet::invalid();
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
        std::array<SizeType, UDim> nn;
        this->dims(nn);
        for(int d = 0; d < Dim; ++d) {
            cell_size[d] = (impl_->box_max_[d] - impl_->box_min_[d])/(nn[d] - 1);
        }
    }

    template<class Elem>
    void FunctionSpace<PetscDM<Elem::Dim>, 1, Elem>::elem(const SizeType &idx, Elem &e) const
    {
        mesh_->elem(idx, e);
        typename Mesh::Point translation, cell_size;
        mesh_->cell_point(idx, translation);
        mesh_->cell_size(idx, cell_size);
        e.set(translation, cell_size);
    }

    template<class Elem>
    bool FunctionSpace<PetscDM<Elem::Dim>, 1, Elem>::write(const Path &path, const PetscVector &x) const
    {
        PetscErrorCode ierr = 0;

        PetscViewer       viewer;
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, path.c_str(), &viewer);
        if(ierr != 0) return false;

        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK); assert(ierr == 0);

        // PetscViewerBinaryOpen(PETSC_COMM_WORLD, path.c_str(),FILE_MODE_WRITE, &viewer);
        // PetscViewerPushFormat(viewer, PETSC_VIEWER_BINARY_VTK);

        DMView(raw_type(*mesh_), viewer);
        VecView(raw_type(x), viewer);

        //Extra output
        PetscVector temp = x;
        temp.set(0.0);
        temp.set(x.comm().rank());

        utopia::rename("comm_rank", temp);
        VecView(raw_type(temp), viewer);


        // each_write(temp, [](const SizeType &i) {
        //     return i;
        // });

        // utopia::rename("global_node_id", temp);
        // VecView(raw_type(temp), viewer);

        ierr = PetscViewerDestroy(&viewer); assert(ierr == 0);
        return ierr == 0;
    }

    template<int Dim>
    void PetscDM<Dim>::node(const SizeType &idx, Point &p) const
    {
        std::array<SizeType, Dim> tensor_index, dims;
        impl_->global_node_grid_coord(idx, tensor_index);

        this->dims(dims);

        for(int i = 0; i < Dim; ++i) {
            SizeType c = tensor_index[i];
            p[i] = c * (impl_->box_max_[i] - impl_->box_min_[i])/(dims[i] - 1) + impl_->box_min_[i];
        }
    }

}

#endif //UTOPIA_PETSC_DM_IMPL_HPP
