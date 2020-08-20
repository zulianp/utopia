
#include "petscts.h"
#include "petsc.h"
int main()
{
  PetscErrorCode ierr;
  TS ts;
  int argc = 0;
  char** argv = NULL;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 1
  ierr = TSDestroy(ts);CHKERRQ(ierr);
#else
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
#endif
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
