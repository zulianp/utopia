// #ifndef UTOPIA_PETSC_FE_HPP
// #define UTOPIA_PETSC_FE_HPP

// #include <vector>

// #include <petscdmplex.h>
// #include <petscdt.h>
// #include <petscfe.h>

// namespace utopia {

//     class PetscFiniteElement {
//     public:
//         using SizeType = PetscInt;

//         void reinit(const SizeType &cell) {
//             PetscQuadrature quad;
//             PetscFEGetQuadrature(fe_, &quad);
//             DMPlexComputeCellGeometryFEM(dmplex_, cell, quad, &points_[0], &J_[0], &J_inv[0], &J_det[0]);
//         }

//         void init(const SizeType &cell) {
//             PetscInt dim;
//             DMGetDimension(dmplex_, &dim);

//             PetscQuadrature quad = nullptr;
//             PetscFEGetQuadrature(fe_, &quad);

//             if (quad) {
//                 SizeType nqp = quad->numPoints;

//                 points_.resize(nqp * dim);
//                 J_.resize(nqp * dim * dim);
//                 J_inv.resize(nqp * dim * dim);
//                 J_det.resize(nqp);

//             } else {
//                 points_.resize(dim);
//                 J_.resize(dim * dim);
//                 J_inv.resize(dim * dim);
//                 J_det.resize(1);
//             }

//             reinit(cell);
//         }

//         PetscFiniteElement(DM dmplex, PetscFE fe) : dmplex_(dmplex), fe_(fe) {}

//     private:
//         DM dmplex_;
//         PetscFE fe_;
//         std::vector<PetscReal> points_;
//         std::vector<PetscReal> J_;
//         std::vector<PetscReal> J_inv;
//         std::vector<PetscReal> J_det;
//     };
// }  // namespace utopia

// #endif  // UTOPIA_PETSC_FE_HPP
