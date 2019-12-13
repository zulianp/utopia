#include "utopia_PetscDM.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_Readable.hpp"

#include <petscdm.h>
#include <petscdmda.h>

#include <csignal>

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

    //https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex42.c.html
    //https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/FE/index.html

    DM raw_type(const PetscDM::Impl &dm);

    class PetscDMImpl {
    public:
        using SizeType = Traits<PetscVector>::SizeType;
        using Scalar   = Traits<PetscVector>::Scalar;

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

        static void box(const DM &dm, Scalar *min, Scalar *max)
        {
            DMDAGetBoundingBox(dm, min, max);
        }

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


       inline static  SizeType dim(const DM &dm)
       {
            SizeType ret;
            DMGetDimension(dm, &ret);
            return ret;
        }


    };

    class PetscDM::Elements {
    public:
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        using Scalar   = utopia::Traits<PetscVector>::Scalar;

        Elements(const PetscDM::Impl &dm_impl)
        : dm(raw_type(dm_impl))
        {
            DMDAGetElements(dm,&ne,&nc,&e);
        }

        template<class Fun>
        void each(Fun fun) const
        {
            ISLocalToGlobalMapping map;
            DMGetLocalToGlobalMapping(dm,&map);
            std::vector<SizeType> e_global(ne*nc);
            ISLocalToGlobalMappingApplyBlock(map, ne*nc, e, &e_global[0]);

            for(SizeType i = 0; i < ne; i++) {
                PetscDM::Elem e(*this, i);
                fun(e);
            }

            PetscSynchronizedFlush(PETSC_COMM_WORLD,stdout);
        }

        ~Elements()
        {
            DMDARestoreElements(dm,&ne,&nc,&e);
        }

        DM dm;
        SizeType ne,nc;
        const SizeType   *e;
    };

    class PetscDM::Nodes {
    public:

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

            each(dim(), dims, start, extent, f);
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

            each(dim(), dims, start, extent, f);
        }

        template<class Fun>
        void each(
            const int dim,
            const SizeType dims[3],
            const SizeType start[3],
            const SizeType extent[3], Fun f) const
        {
            switch(dim) {
                case 1:
                {
                    for(SizeType i = 0; i < extent[0]; ++i) {
                        SizeType idx = i + start[0];

                        PetscDM::Node node(*this, idx);
                        f(node);
                    }

                    break;
                }

                case 2:
                {
                    for(SizeType j = 0; j < extent[1]; ++j) {
                        for(SizeType i = 0; i < extent[0]; ++i) {
                            SizeType idx = (j + start[1]) * dims[0] + i + start[0];

                            PetscDM::Node node(*this, idx);
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

                                PetscDM::Node node(*this, idx);
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

        Nodes(const PetscDM::Impl &impl)
        : dm(raw_type(impl)), range(PetscDMImpl::local_node_range(dm))
        {}

    private:
        DM dm;
        Range range;
    };

    class PetscDM::Impl {
    public:
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        using Scalar   = utopia::Traits<PetscVector>::Scalar;

        Impl()
        : dm(nullptr)
        {}

        inline Range local_node_range() const
        {
            return PetscDMImpl::local_node_range(dm);
        }

        void create_matrix(PetscMatrix &mat)
        {
            mat.destroy();
            DMCreateMatrix(dm, &mat.raw_type());
        }

        void create_vector(PetscVector &vec)
        {
            vec.destroy();
            DMCreateGlobalVector(dm, &vec.raw_type());
        }

        template<class Fun>
        void each_element(Fun fun) const
        {
            Elements elements(*this);
            elements.each(fun);
        }

        template<class Fun>
        void each_node(Fun fun) const
        {
            Nodes nodes(*this);
            nodes.each(fun);
        }

        template<class Fun>
        void each_node_with_ghosts(Fun fun) const
        {
            Nodes nodes(*this);
            nodes.each_with_ghosts(fun);
        }

        inline SizeType elem_node_id(const SizeType &e, const SizeType &k) const
        {
            return -1;
        }

        inline SizeType elem_n_nodes() const
        {
            return std::pow(dim(), 2);
        }

        template<class Array>
        void init(
            const PetscCommunicator &comm,
            const Array &arr,
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

        template<class Array>
        void dims(Array &arr) const
        {
            auto ierr = DMDAGetInfo(dm, nullptr, &arr[0], &arr[1], &arr[2], nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); assert(ierr == 0);
        }

        template<class Array, class Coords>
        void init_uniform(
            const PetscCommunicator &comm,
            const Array &arr,
            const Coords &box_min,
            const Coords &box_max,
            DMDAElementType elem_type    = DMDA_ELEMENT_Q1,
            DMDAStencilType stencil_type = DMDA_STENCIL_BOX)
        {
            init(comm, arr, stencil_type);

            Scalar min_x = 0, min_y = 0, min_z = 0;
            Scalar max_x = 0, max_y = 0, max_z = 0;

            min_x = box_min[0];
            max_x = box_max[0];

            const auto d = arr.size();
            if(d > 1) {
                min_y = box_min[1];
                max_y = box_max[1];
            }

            if(d > 2) {
                min_z = box_min[2];
                max_z = box_max[2];
            }

            DMDASetUniformCoordinates(dm, min_x, max_x, min_y, max_y, min_z, max_z);
            DMDASetElementType(dm, elem_type);
            DMSetUp(dm);
        }

        ~Impl() {
            destroy();
        }

        void destroy()
        {
            if(dm) {
                DMDestroy(&dm);
                dm = nullptr;
            }
        }

        void describe() const
        {

            PetscErrorCode err = 0;

            // DMDALocalInfo info;
            // err = DMDAGetLocalInfo(dm,&info); assert(err == 0);
            // PetscViewer viewer;
            // PetscViewerDrawOpen(comm.get(),0,"",300,0,300,300,&viewer);
            // DMView(dm, viewer);
            // std::raise(SIGINT);
            // PetscViewerDestroy(&viewer);

            SizeType dim, M, N, P, m, n, p, dof, s;
            DMDAStencilType st;
            err = DMDAGetInfo(dm, &dim, &M, &N, &P, &m, &n, &p, &dof, &s, nullptr, nullptr, nullptr, &st); assert(err == 0);

            int size = comm.size();
            int rank = comm.rank();

            comm.barrier();

            Vec v;
            DMGetGlobalVector(dm, &v);

            PetscInt n_b, n_e;
            VecGetOwnershipRange(v, &n_b, &n_e);
            Range range(n_b, n_e);

            for(int i = 0; i < size; ++i) {
                if(i == rank) {
                    std::cout << "--------------------------------------------\n";
                    std::cout << i << ")\ndim: " << dim << "\nM: " << M << "\nN: " << N << "\nP: " << P << "\nm: " << m << "\nn: " << n << "\np: " << p << "\ndof: " << dof << "\ns: " << s <<  std::endl;

                    Nodes nodes(*this);

                    nodes.each_with_ghosts([range](const PetscDM::Node &node) {

                        std::cout << node.idx() << " ";
                        if(range.inside(node.idx())) {

                        } else {
                            std::cout << "g";
                        }

                        std::cout << std::endl;
                    });


                    std::cout << "--------------------------------------------\n";
                }

                comm.barrier();
            }

            DMRestoreGlobalVector(dm, &v);

            // print_ghost_index();


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


        // void print_ghost_index() const
        // {
        //     IS is;
        //     DMDAGetSubdomainCornersIS(dm, &is);

        //     PetscIS idx(is, [this](IS *is) -> PetscErrorCode {
        //        return DMDARestoreSubdomainCornersIS(dm, is);
        //     });


        //     ISLocalToGlobalMapping map;
        //     DMGetLocalToGlobalMapping(dm,&map);

        //     Read<PetscIS> ri(idx);

        //     int size = comm.size();
        //     int rank = comm.rank();

        //     comm.barrier();

        //     PetscInt M, N, P;
        //     PetscInt p_sizes[3] = {0, 0, 0};
        //     const PetscInt *ranges[3];
        //     DMDAGetOwnershipRanges(dm, &ranges[0], &ranges[1], &ranges[2]);

        //     DMDAGetInfo(
        //         dm, 0,
        //         &M, &N, &P,
        //         &p_sizes[0], &p_sizes[1], &p_sizes[3],
        //         0,
        //         0,
        //         0,0,0,
        //         0
        //     );


        //     Vec v;
        //     DMGetGlobalVector(dm, &v);

        //     PetscInt n_b, n_e;
        //     VecGetOwnershipRange(v, &n_b, &n_e);

        //     for(int i = 0; i < size; ++i) {
        //         if(i == rank) {
        //             std::cout << "--------------------------------------------\n";

        //             auto n = idx.size();
        //             std::vector<SizeType> gidx(n);

        //             ISLocalToGlobalMappingApply(map, n, &idx[0], &gidx[0]);
        //             for(SizeType k = 0; k < n; ++k) {
        //                 std::cout << gidx[k] << " ";
        //             }

        //             std::cout << std::endl;
        //             std::cout << n_b << " " << n_e << std::endl;
        //             std::cout << std::endl;
        //             std::cout << "--------------------------------------------\n";
        //         }

        //         comm.barrier();
        //     }

        //     DMRestoreGlobalVector(dm, &v);
        // }

        inline SizeType dim() const
        {
            SizeType ret;
            DMGetDimension(dm, &ret);
            return ret;
        }

        DM dm;
        PetscCommunicator comm;
    };

    PetscDM::PetscDM()
    : impl_(utopia::make_unique<Impl>())
    {}

    DM raw_type(const PetscDM::Impl &dm)
    {
        return dm.dm;
    }

    Range local_node_range(const PetscDM::Impl &dm)
    {
        return dm.local_node_range();
    }

    PetscDM::PetscDM(
        const PetscCommunicator &comm,
        const std::vector<SizeType> &dims,
        const std::vector<Scalar> &box_min,
        const std::vector<Scalar> &box_max)
    : impl_(utopia::make_unique<Impl>())
    {
        impl_->init_uniform(comm, dims, box_min, box_max);
    }

    void PetscDM::describe() const
    {
        impl_->describe();
    }

    PetscDM::~PetscDM()
    {

    }

    void PetscDM::create_matrix(PetscMatrix &mat)
    {
        impl_->create_matrix(mat);
    }

    PetscDM::SizeType PetscDM::Elem::n_nodes() const
    {
        // return impl_.elem_n_nodes();
        return -1;
    }

    PetscDM::SizeType PetscDM::Elem::node_id(const SizeType k) const
    {
        // return impl_.elem_node_id(idx_, k);
        return -1;
    }

    void PetscDM::each_element(const std::function<void(const Elem &)> &f)
    {
        impl_->each_element(f);
    }

    void PetscDM::each_node(const std::function<void(const Node &)> &f)
    {
        impl_->each_node(f);
    }

    void PetscDM::each_node_with_ghosts(const std::function<void(const Node &)> &f)
    {
        impl_->each_node_with_ghosts(f);
    }

    bool PetscDM::Node::is_ghost() const
    {
        return nodes_.is_ghost(idx());
    }

    Range PetscDM::local_node_range() const
    {
        return impl_->local_node_range();
    }

    void PetscDM::dims(SizeType *arr) const
    {
       impl_->dims(arr);
    }

    PetscDM::SizeType PetscDM::dim() const
    {
        return impl_->dim();
    }

    void PetscDM::box(Scalar *min, Scalar *max) const
    {
        PetscDMImpl::box(impl_->dm, min, max);
    }

    void PetscDM::local_node_ranges(SizeType *begin, SizeType *end) const
    {
        PetscDMImpl::local_node_ranges(impl_->dm, begin, end);
    }

    void PetscDM::local_element_ranges(SizeType *begin, SizeType *end) const
    {
        PetscDMImpl::local_element_ranges(impl_->dm, begin, end);
    }

    PetscDM::SizeType PetscDM::n_local_nodes_with_ghosts() const
    {
        return PetscDMImpl::n_local_nodes_with_ghosts(impl_->dm);
    }
}
