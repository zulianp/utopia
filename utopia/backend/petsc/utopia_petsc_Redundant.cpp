#include "utopia_petsc_Redundant.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    class PetscVecScatter::Wrapper {
    public:
        Wrapper() : scatter(nullptr), owned_(true) {}

        void destroy()
        {
            if(owned_ && scatter)
            {
                VecScatterDestroy(&scatter);
                scatter = nullptr;
                owned_  = true;
            }
        }

        virtual ~Wrapper()
        {
            destroy();
        }

        inline VecScatter &raw_type()
        {
            return scatter;
        }

    private:
        VecScatter scatter;
        bool owned_;
    };

    PetscVecScatter::PetscVecScatter() {}
    PetscVecScatter::~PetscVecScatter() {}

    void PetscVecScatter::create(
        const PetscVector &from, const PetscIS &from_is,
        const PetscVector &to, const PetscIS &to_is
    )
    {
        if(!wrapper_) {
            wrapper_ = std::make_shared<Wrapper>();
        } else {
            wrapper_->destroy();
        }

        VecScatterCreate(from.raw_type(), from_is.raw_type(), to.raw_type(), to_is.raw_type(), &wrapper_->raw_type());
    }

    void PetscVecScatter::apply(const PetscVector &from, PetscVector &to) const
    {
        begin(from, to);
        end(from, to);
    }

    void PetscVecScatter::begin(const PetscVector &from, PetscVector &to) const
    {
        VecScatterBegin(wrapper_->raw_type(), from.raw_type(), to.raw_type(), INSERT_VALUES, SCATTER_FORWARD);
    }

    void PetscVecScatter::end(const PetscVector &from, PetscVector &to) const
    {
        VecScatterEnd(wrapper_->raw_type(), from.raw_type(), to.raw_type(), INSERT_VALUES, SCATTER_FORWARD);
    }

    Redundant<PetscMatrix, PetscVector>::Redundant()
    : n_sub_comm_(2)
    {}

    Redundant<PetscMatrix, PetscVector>::~Redundant()
    {
        // VecScatterDestroy(&scatterout);
        // VecScatterDestroy(&scatterin);
        PetscSubcommDestroy(&psubcomm);
    }

    // Redundant<PetscMatrix, PetscVector> * Redundant<PetscMatrix, PetscVector>::clone() const
    // {
    //     return new Redundant(*this);
    // }

    // construct idx sets and paritions
    void Redundant<PetscMatrix, PetscVector>::init(const Layout &lo, const SizeType n_sub_comm)
    {
        n_sub_comm_ = n_sub_comm;

        assert(n_sub_comm > 1);

        auto comm = lo.comm().get();

        //create sub-comm
        PetscSubcommCreate(comm, &psubcomm);
        PetscSubcommSetNumber(psubcomm, n_sub_comm_);
        PetscSubcommSetType(psubcomm, PETSC_SUBCOMM_CONTIGUOUS);
        MPI_Comm subcomm = PetscSubcommChild(psubcomm);




        //create paritioning
        PetscInt start_sub = 0, end_sub = 0, local_size_sub = PETSC_DECIDE, size = lo.size();
        PetscSplitOwnership(subcomm, &local_size_sub, &size);
        MPI_Scan(&local_size_sub, &end_sub, 1, MPIU_INT, MPI_SUM, subcomm);
        start_sub = end_sub - local_size_sub;

        PetscInt size_sub = n_sub_comm * size;

        sub_layout_ = layout(PetscCommunicator(PetscSubcommContiguousParent(psubcomm)), local_size_sub, size_sub);

        int mpi_compare;
        MPI_Comm_compare(comm, sub_layout_.comm().get(), &mpi_compare);

        assert(mpi_compare != MPI_UNEQUAL); //MPI_SIMILAR

        // create index-sets
        PetscInt *idx1, *idx2;

        PetscMalloc2(psubcomm->n*local_size_sub, &idx1, psubcomm->n*local_size_sub, &idx2);

        PetscInt j = 0;

        for(PetscInt k=0; k < psubcomm->n; k++) {
            for (PetscInt i = start_sub; i < end_sub; i++) {
                idx1[j]   = i;
                idx2[j++] = i + size * k;
            }
        }

        {
            IS is1, is2;
            ISCreateGeneral(comm, psubcomm->n*local_size_sub, idx1, PETSC_OWN_POINTER, &is1);
            ISCreateGeneral(comm, psubcomm->n*local_size_sub, idx2, PETSC_OWN_POINTER, &is2);
            is_super_to_sub_from = utopia::make_unique<PetscIS>(is1);
            is_super_to_sub_to   = utopia::make_unique<PetscIS>(is2);
        }

        {
            IS is1, is2;
            ISCreateStride(comm, local_size_sub, start_sub + psubcomm->color * size, 1, &is1);
            ISCreateStride(comm, local_size_sub, start_sub, 1, &is2);
            is_sub_to_super_from = utopia::make_unique<PetscIS>(is1);
            is_sub_to_super_to   = utopia::make_unique<PetscIS>(is2);
        }
    }

    void Redundant<PetscMatrix, PetscVector>::create_sub_vector(
        const PetscVector &vec,
        PetscVector &vec_sub,
        PetscVecScatter &scatter_to_sub,
        PetscVecScatter &scatter_to_super)
    {
        //create vec_sub
        vec_sub.zeros(sub_layout_);

        //create scatters
        scatter_to_sub.create(vec,       *is_super_to_sub_from, vec_sub, *is_super_to_sub_to);
        scatter_to_super.create(vec_sub, *is_sub_to_super_from, vec,     *is_sub_to_super_to);
    }

    void Redundant<PetscMatrix, PetscVector>::super_to_sub(const PetscVector &vec, const PetscVecScatter &scatter, PetscVector &vec_sub)
    {
        scatter.apply(vec, vec_sub);
    }

    void Redundant<PetscMatrix, PetscVector>::sub_to_super(const PetscVector &vec_sub, const PetscVecScatter &scatter, PetscVector &vec)
    {
        scatter.apply(vec_sub, vec);
    }

    void Redundant<PetscMatrix, PetscVector>::create_sub_matrix(const PetscMatrix &mat, PetscMatrix &mat_sub)
    {
        MPI_Comm subcomm = PetscSubcommChild(psubcomm);

        mat_sub.destroy();
        mat_sub.comm().set(subcomm);

        MatCreateRedundantMatrix(
            mat.raw_type(),
            psubcomm->n,
            subcomm,
            MAT_INITIAL_MATRIX,
            &mat_sub.raw_type()
        );
    }

    void Redundant<PetscMatrix, PetscVector>::super_to_sub(const PetscMatrix &mat, PetscMatrix &mat_sub)
    {
        MatCreateRedundantMatrix(mat.raw_type(), psubcomm->n, PetscSubcommChild(psubcomm), MAT_REUSE_MATRIX, &mat_sub.raw_type());
    }

}
