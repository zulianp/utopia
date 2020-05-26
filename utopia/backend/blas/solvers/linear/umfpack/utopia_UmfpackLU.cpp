// FIXME
// #include "utopia_UmfpackLU.hpp"
// #include <umfpack.h>

// namespace utopia {

// namespace internals {
//     UmfpackLU::~UmfpackLU()
//     {

//     }

//     bool UmfpackLU::solve(const CCSMatrixd &mat, const Vectord &rhs, Vectord &solution)
//     {
//         if(solution.size() != rhs.size()) {
//             solution = zeros(rhs.size());
//         }

//         Internal internal;
//         init(mat, internal);
//         bool ok = solve(mat, internal, ptr(rhs), ptr(solution));
//         cleanUp(internal);
//         return ok;
//     }

//     bool UmfpackLU::solve(const CRSMatrixd &matrix, const Vectord &rhs, Vectord &solution)
//     {
//         CCSMatrixd ccsmat = matrix;
//         return solve(ccsmat, rhs, solution);
//     }

//     void UmfpackLU::cleanUp(Internal &internal) const
//     {
//         if(internal.numeric) {
//             umfpack_di_free_numeric(&internal.numeric);
//         }
//     }

//     bool UmfpackLU::init(const CCSMatrixd &mat, Internal &internal) const
//     {
//         cleanUp(internal);

//         std::vector<double> &info = internal.info;
//         std::vector<double> &control = internal.control;

//         int status = 0;

//         control.resize(UMFPACK_CONTROL);
//         info.resize(UMFPACK_INFO);

//         umfpack_di_defaults(&control[0]);

//         status = umfpack_di_symbolic(mat.size().get(0),                          //0
//                                      mat.size().get(1),                          //1
//                                      &mat.implementation().colptr()[0],          //2
//                                      &mat.implementation().rowindex()[0],        //3
//                                      &mat.implementation().entries()[0],      //4
//                                      &internal.symbolic, &control[0], &info[0]);       //5, 6, 7

//         assert(status == UMFPACK_OK);

//         status = umfpack_di_numeric (&mat.implementation().colptr()[0],          //0
//                                      &mat.implementation().rowindex()[0],        //1
//                                      &mat.implementation().entries()[0],      //2
//                                      internal.symbolic, &internal.numeric, &control[0], &info[0]);     //4, 5, 6, 7

//         assert(status == UMFPACK_OK);

//         if(internal.symbolic) {
//             umfpack_di_free_symbolic(&internal.symbolic);
//         }

//         return (status == UMFPACK_OK);
//     }

//     bool UmfpackLU::solve(const CCSMatrixd &mat, Internal &internal, const double * rhs, double * solution)
//     {
//         int status = umfpack_di_solve (UMFPACK_A,
//                                        &mat.implementation().colptr()[0],
//                                        &mat.implementation().rowindex()[0],
//                                        &mat.implementation().entries()[0],
//                                        solution, rhs, internal.numeric,
//                                        &internal.control[0],
//                                        &internal.info[0]);

//         if(status != UMFPACK_OK) {
//             std::cerr << "status != UMFPACK_OK --> " << status << std::endl;
//             disp(mat, std::cerr);
//         }
//         assert(status == UMFPACK_OK);

//         return (status == UMFPACK_OK);
//     }

// }
// }