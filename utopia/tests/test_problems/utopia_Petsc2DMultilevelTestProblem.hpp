#include "utopia.hpp"

#ifdef  WITH_PETSC
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscmatlab.h>
#include <petsc/private/snesimpl.h> /* For SNES_Solve event */

namespace utopia
{
    template<typename Matrix, typename Vector, typename ProblemType>
    class Petsc2DMultilevelTestProblem
    {
        public:
            typedef UTOPIA_SIZE_TYPE(DVectord) SizeType;
            typedef UTOPIA_SCALAR(DVectord) Scalar;


        Petsc2DMultilevelTestProblem(const SizeType &n_levels = 2, const SizeType & n_coarse = 10, bool remove_BC_contributions = false):
            n_levels_(n_levels),
            n_coarse_(n_coarse),
            remove_BC_contributions_(remove_BC_contributions)
        {
            std::vector<DM> dms_; 
            dms_.resize(n_levels); 
            level_functions_.resize(n_levels); 

            DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR, n_coarse_, n_coarse_, PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL, &dms_[0]);
            DMSetUp(dms_[0]);
            DMDASetUniformCoordinates(dms_[0], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

            for(auto l=1; l <n_levels; l++)
            {
                DMRefine(dms_[l-1],PETSC_COMM_WORLD, &dms_[l]); 
            }

            for(auto l=1; l <n_levels; l++)
            {
                Mat I;
                Matrix I_u; 
                DMCreateInterpolation(dms_[l-1], dms_[l], &I, 0);
                wrap(I, I_u);
                transfers_.push_back( std::make_shared<IPTransfer<Matrix, Vector> >( std::make_shared<Matrix>(I_u), 0.5) );
                MatDestroy(&I);
            }

            for(auto l=0; l <n_levels; l++)
            {
                auto fun = std::make_shared<ProblemType >(dms_[l]);
                level_functions_[l] = fun;           
            }                

        }     

        ~Petsc2DMultilevelTestProblem()
        {
            level_functions_.clear(); 
            transfers_.clear(); 
        }


        SizeType n_levels_;
        SizeType n_coarse_;
        bool remove_BC_contributions_;

        std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > transfers_;
        std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions_; 

    };
}

#endif //WITH_PETSC
