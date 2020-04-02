#include "utopia_petsc_Redundant.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_PetscIS.hpp"

namespace utopia {

    class PetscVecScatter {
    public:



    private:
        class Wrapper {
        public:
            Wrapper() : scatter(nullptr), owned_(true) {}

            virtual ~Wrapper()
            {
                if(owned_ && scatter)
                {
                    VecScatterDestroy(&scatter);
                }
            }

            VecScatter scatter;
            bool owned_;
        };

        std::shared_ptr<Wrapper> wrapper_;
    };

    Redundant<PetscMatrix, PetscVector>::Redundant()
    : n_sub_comm_(2)
    {}

    Redundant<PetscMatrix, PetscVector>::~Redundant()
    {
        VecScatterDestroy(&scatterout);
        VecScatterDestroy(&scatterin);
        PetscSubcommDestroy(&psubcomm);
    }

    Redundant<PetscMatrix, PetscVector> * Redundant<PetscMatrix, PetscVector>::clone() const
    {
        return new Redundant(*this);
    }

    // construct idx sets and paritions
    void Redundant<PetscMatrix, PetscVector>::init(const Layout &lo)
    {
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

        // create index-sets
        IS        is1, is2;
        PetscInt *idx1, *idx2;

        PetscMalloc2(psubcomm->n*local_size_sub, &idx1, psubcomm->n*local_size_sub, &idx2);

        PetscInt j = 0;

        for(PetscInt k=0; k < psubcomm->n; k++) {
            for (PetscInt i=start_sub; i< end_sub; i++) {
                idx1[j]   = i;
                idx2[j++] = i + size * k;
            }
        }

        // ISCreateGeneral(comm, psubcomm->n*local_size_sub, idx1, PETSC_COPY_VALUES, &is1);
        // ISCreateGeneral(comm, psubcomm->n*local_size_sub, idx2, PETSC_COPY_VALUES, &is2);
        // VecScatterCreate(raw_type(x), is1, raw_type(sol_dup), is2, &scatterin);
        // ISDestroy(&is1);
        // ISDestroy(&is2);

        // ISCreateStride(comm, local_size_sub, mstart+ psubcomm->color*M, 1, &is1);
        // ISCreateStride(comm, local_size_sub, mstart, 1, &is2);
        // VecScatterCreate(raw_type(sol_dup), is1, raw_type(x), is2, &scatterout);
        // ISDestroy(&is1);
        // ISDestroy(&is2);
        PetscFree2(idx1, idx2);
    }

    void Redundant<PetscMatrix, PetscVector>::create_sub_vector(
        const PetscVector &vec,
        PetscVector &vec_sub,
        PetscVecScatter &scatter_to_sub,
        PetscVecScatter &scatter_to_super)
    {

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

    // void Redundant<PetscMatrix, PetscVector>::super_to_sub(const PetscMatrix &mat, const PetscVecScatter &scatter, PetscMatrix &mat_sub)
    // {

    // }

    // void Redundant<PetscMatrix, PetscVector>::super_to_sub(const PetscVector &vec, const PetscVecScatter &scatter, PetscVector &vec_sub)
    // {
    // // scatter vecution, rhs
    //     // PetscScalar    *array_vec, *array_rhs;

    //     // VecScatterBegin(scatterin, raw_type(vec), raw_type(vec_dup), INSERT_VALUES, SCATTER_FORWARD);
    //     // VecScatterEnd(scatterin, raw_type(vec), raw_type(vec_dup), INSERT_VALUES, SCATTER_FORWARD);

    //     // VecGetArray(raw_type(vec_dup), &array_vec);
    //     // VecPlaceArray(raw_type(vec_sub),  (const PetscScalar*)array_vec);
    // }

    void Redundant<PetscMatrix, PetscVector>::sub_to_super(const PetscVector &vec_sub, PetscVector &vec)
    {
        // PetscScalar    *array_vec;

        // /* place ysub's local array into ydup */
        // VecGetArray(raw_type(vec_sub), &array_vec);
        // VecPlaceArray(raw_type(vec_dup), (const PetscScalar*)array_vec);


        // /* scatter ydup to y */
        // VecScatterBegin(scatterout, raw_type(vec_dup), raw_type(vec), INSERT_VALUES, SCATTER_FORWARD);
        // VecScatterEnd(scatterout, raw_type(vec_dup), raw_type(vec), INSERT_VALUES, SCATTER_FORWARD);


        // VecResetArray(raw_type(vec_dup));
        // VecRestoreArray(raw_type(vec_sub), &array_vec);
    }

}
