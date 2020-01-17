#ifndef UTOPIA_LINEAR_SOLVE_EXPORT_HPP
#define UTOPIA_LINEAR_SOLVE_EXPORT_HPP

typedef struct {
    void *ptr;
} USolverImpl;

typedef struct {
    void *ptr;
} UMatImpl;

typedef struct {
    void *ptr;
} UVecImpl;

typedef USolverImpl * USolver;
typedef UMatImpl    * UMat;
typedef UVecImpl    * UVec;

typedef const char * USolverType;
typedef const char * UPreconditionerType;
typedef const char * UPackage;

static USolverType U_KSP = "ksp";
static UPackage    U_BJACOBI = "bjacobi";
static USolverType U_PETSC = "petsc";

void UtopiaInitialize(int argc, char *argv[]);
int UtopiaFinalize();
void USolverCreate(USolver * solver, USolverType type, UPreconditionerType prec, UPackage package);
void USolverDestroy(USolver * solver);
void USolverSolve(USolver solver, UMat A, UVec b, UVec x);
void USolverPrintInfo(USolver ptr);
void UtopiaPrintVersion();

#endif
