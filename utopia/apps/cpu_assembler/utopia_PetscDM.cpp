#include "utopia_PetscDM.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Communicator.hpp"

#include <petscdm.h>
#include <petscdmda.h>

#include <csignal>

namespace utopia {

    //https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex42.c.html

    DM raw_type(const PetscDM::Impl &dm);

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
            // PetscDM::Elem e(*this, 0);

            ISLocalToGlobalMapping map;
            DMGetLocalToGlobalMapping(dm,&map);
            std::vector<SizeType> e_global(ne*nc);
            ISLocalToGlobalMappingApplyBlock(map, ne*nc, e, &e_global[0]);

            for(SizeType i = 0; i < ne; i++) {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"i %D e[4*i] %D %D %D %D\n",i,
                    // e[4*i],e[4*i+1],e[4*i+2], e[4*i+3]
                    e_global[4*i],e_global[4*i+1],e_global[4*i+2], e_global[4*i+3]
                );
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

    class PetscDM::Impl {
    public:
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        using Scalar   = utopia::Traits<PetscVector>::Scalar;

        Impl()
        : dm(nullptr)
        {}

        void create_matrix(PetscMatrix &mat)
        {
            mat.destroy();
            DMCreateMatrix(dm, &mat.raw_type());
        }



        template<class Fun>
        void each_element(Fun fun) const
        {
            PetscErrorCode ierr;
            Elements elements(*this);
            elements.each(fun);
        }

        inline SizeType elem_node_id(const SizeType &e, const SizeType &k) const
        {
            return -1;
        }

        inline SizeType elem_n_nodes() const
        {
            return std::pow(dim, 2);
        }

        template<class Array>
        void init(
            const PetscCommunicator &comm,
            const Array &arr)
        {
            destroy();
            this->comm = comm;

            dim = arr.size();

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
                        DMDA_STENCIL_BOX,
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
                        DMDA_STENCIL_BOX,
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

        template<class Array, class Coords>
        void init_uniform(
            const PetscCommunicator &comm,
            const Array &arr,
            const Coords &box_min,
            const Coords &box_max)
        {

            init(comm, arr);


            Scalar min_x = 0, min_y = 0, min_z = 0;
            Scalar max_x = 0, max_y = 0, max_z = 0;

            min_x = box_min[0];
            max_x = box_max[0];

            if(dim > 1) {
                min_y = box_min[1];
                max_y = box_max[1];
            }

            if(dim > 2) {
                min_z = box_min[2];
                max_z = box_max[2];
            }

            DMDASetUniformCoordinates(dm, min_x, max_x, min_y, max_y, min_z, max_z);
            // DMDASetElementType(dm, DMDA_ELEMENT_P1);
            DMDASetElementType(dm, DMDA_ELEMENT_Q1);
            DMSetUp(dm);

            init_aux();
        }

        ~Impl() {
            destroy();
        }

        void destroy()
        {
            if(dm) {
                DMDestroy(&dm);
                dm = nullptr;
                dim = 0;
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

            for(int i = 0; i < size; ++i) {
                if(i == rank) {
                    std::cout << "--------------------------------------------\n";
                    std::cout << i << ")\ndim: " << dim << "\nM: " << M << "\nN: " << N << "\nP: " << P << "\nm: " << m << "\nn: " << n << "\np: " << p << "\ndof: " << dof << "\ns: " << s <<  std::endl;

                    for(SizeType d = 0; d < dim; ++d) {
                        std::cout << "[" << corners_start[d] << ", " << corners_end[d] << "]" << " g = ";
                        std::cout << "[" << ghost_corners_start[d] << ", " << ghost_corners_end[d] << "]" << std::endl;
                    }

                    std::cout << "--------------------------------------------\n";
                }

                comm.barrier();
            }

            print_ghost_index();

        }


        void print_ghost_index() const
        {
            ISLocalToGlobalMapping map;
            DMGetLocalToGlobalMapping(dm,&map);

            std::vector<SizeType> coord(dim);
            for(SizeType d = 0; d < dim; ++d) {
                coord[d] = corners_end[d];

                const SizeType end_layers   = ghost_corners_end[d] - corners_end[d];
                const SizeType begin_layers = corners_end[d] - ghost_corners_begin[d];

                std::cout << "layers(" << d << ") " <<begin_layers + end_layers << std::endl;

                if(end_layers == 0) continue;

                std::cout <<

                for(SizeType d2 = 0; d2 < dim; ++d2) {
                    if(d2 == d) continue;

                    if(dim == 2) {
                        const SizeType n2 = ghost_corners_end[d2] - corners_end[d2];
                        for(SizeType j = 0; j < n2; ++j) {
                            coord[d2] = j;

                            std::cout << "(" << coord[0] << " " << coord[1] << ")" << std::endl;
                        }
                    }

                }

            }

            // ISLocalToGlobalMappingApply(map, PetscInt N,const PetscInt in[],PetscInt out[])
        }

        DM dm;
        SizeType dim;
        PetscCommunicator comm;

        //array version of dm's data
        std::vector<SizeType> corners_start;
        std::vector<SizeType> corners_end;

        std::vector<SizeType> ghost_corners_start;
        std::vector<SizeType> ghost_corners_end;

        std::vector<SizeType> local_to_global;

    private:
        void init_aux() {
            corners_start.resize(3, 0);
            corners_end.resize(3, 0);
            ghost_corners_start.resize(3, 0);
            ghost_corners_end.resize(3, 0);

            DMDAGetCorners(dm,
                &corners_start[0], &corners_start[1], &corners_start[2],
                &corners_end[0],   &corners_end[1],   &corners_end[2]
            );

            DMDAGetGhostCorners(dm,
                &ghost_corners_start[0], &ghost_corners_start[1], &ghost_corners_start[2],
                &ghost_corners_end[0],   &ghost_corners_end[1],   &ghost_corners_end[2]
            );

            corners_start.resize(dim);
            corners_end.resize(dim);
            ghost_corners_start.resize(dim);
            ghost_corners_end.resize(dim);
        }
    };


    PetscDM::PetscDM()
    : impl_(utopia::make_unique<Impl>())
    {}

    DM raw_type(const PetscDM::Impl &dm)
    {
        return dm.dm;
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

    void PetscDM::each_element(const std::function<void(const SizeType &, const Elem &e)> &f)
    {
        impl_->each_element(f);
    }

}
