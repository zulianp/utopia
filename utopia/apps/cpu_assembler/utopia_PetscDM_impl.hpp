#ifndef UTOPIA_PETSC_DM_IMPL_HPP
#define UTOPIA_PETSC_DM_IMPL_HPP


#include "utopia_PetscDM.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_Readable.hpp"
#include "utopia_PetscIS.hpp"

#include "utopia_UniformQuad4.hpp"

#include <petscdm.h>
#include <petscdmda.h>

#include <csignal>

namespace utopia {

    class PetscUniformQuad4 final : public PetscElem<2> {
    public:
        using Super    = utopia::PetscElem<2>;
        using SizeType = Super::SizeType;
        using Scalar   = Super::Scalar;
        using Point    = Super::Point;
        using Grad     = Super::Grad;

        inline Scalar fun(const SizeType &i, const Point &p) const
        {
            return UniformQuad4<Scalar>::fun(i, p);
        }

        inline void grad(const int i, const Point &p, Grad &g) const
        {
            RefQuad4::grad(i, p, g);
            g(0) /= h_(0);
            g(1) /= h_(1);
        }

        inline constexpr static bool is_affine()
        {
            return true;
        }

        inline constexpr static Scalar reference_measure()
        {
            return 1.0;
        }

        inline constexpr static int n_nodes()
        {
            return 4;
        }

        inline void set(
            const StaticVector2<Scalar> &translation,
            const StaticVector2<Scalar> &h)
        {
            translation_(0) = translation(0);
            translation_(1) = translation(1);

            h_(0) = h(0);
            h_(1) = h(1);
        }

    private:
        StaticVector2<Scalar> h_;
        StaticVector2<Scalar> translation_;
    };


    //https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex42.c.html
    //https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/FE/index.html

    template<int Dim>
    DM raw_type(const PetscDM<Dim> &dm);

    template<int Dim>
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
        }

        template<class Fun>
        void each(Fun fun) const
        {
            ISLocalToGlobalMapping map;
            DMGetLocalToGlobalMapping(dm,&map);
            std::vector<SizeType> e_global(ne*nc);
            ISLocalToGlobalMappingApplyBlock(map, ne*nc, e, &e_global[0]);

            PetscUniformQuad4 quad4;

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
    };

    template<int Dim>
    class PetscDMNodes {
    public:
        using SizeType = PetscInt;
        using Scalar   = PetscScalar;
        using PetscNode = utopia::PetscNode<Dim>;

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

        void destroy()
        {
            if(dm) {
                DMDestroy(&dm);
                dm = nullptr;
            }
        }

        PetscCommunicator comm;
        DM dm;
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
    }

    template<int Dim>
    PetscDM<Dim>::~PetscDM()
    {}

    template<int Dim>
    void PetscDM<Dim>::each_element(const std::function<void(const Elem &)> &f)
    {
        PetscDMElements<Dim> elements(*this);
        elements.each(f);
    }

    template<int Dim>
    void PetscDM<Dim>::each_node(const std::function<void(const Node &)> &f)
    {
        PetscDMNodes<Dim> nodes(*this);
        nodes.each(f);
    }

    template<int Dim>
    void PetscDM<Dim>::each_node_with_ghosts(const std::function<void (const Node &)> &f)
    {
        PetscDMNodes<Dim> nodes(*this);
        nodes.each_with_ghosts(f);
    }

    template<int Dim>
    Range PetscDM<Dim>::local_node_range() const
    {
        return PetscDMImpl<Dim>::local_node_range(impl_->dm);
    }

    template<int Dim>
    void PetscDM<Dim>::box(Scalar *min, Scalar *max) const
    {
        PetscDMImpl<Dim>::box(impl_->dm, min, max);
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

    template<int Dim>
    bool PetscNode<Dim>::is_ghost() const
    {
        return nodes_.is_ghost(idx());
    }
}

#endif //UTOPIA_PETSC_DM_IMPL_HPP
