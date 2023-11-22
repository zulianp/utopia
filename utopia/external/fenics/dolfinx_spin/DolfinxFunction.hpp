#ifndef UTOPIA_DOLFINX_FUNCTION_HPP
#define UTOPIA_DOLFINX_FUNCTION_HPP

#include <basix/finite-element.h>
#include <dolfinx.h>
#include <dolfinx/common/log.h>
#include <dolfinx/fem/assembler.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/la/Vector.h>
#include <dolfinx/mesh/Mesh.h>
#include <dolfinx/mesh/cell_types.h>
#include <dolfinx/nls/NewtonSolver.h>
#include <cmath>

#include "utopia.hpp"

class DolfinxFunction : public utopia::Function<utopia::PetscMatrix, utopia::PetscVector> {
public:
	using T = PetscScalar;

    class Impl;

    DolfinxFunction(
    	std::shared_ptr<dolfinx::fem::FunctionSpace> V,
             std::shared_ptr<dolfinx::fem::Function<T>> u,
             std::shared_ptr<dolfinx::fem::Form<T>> objective,
             std::shared_ptr<dolfinx::fem::Form<T>> gradient,
             std::shared_ptr<dolfinx::fem::Form<T>> hessian,
             std::vector<std::shared_ptr<const dolfinx::fem::DirichletBC<T>>> boundary_conditions);
    virtual ~DolfinxFunction();

    void create_vector(utopia::PetscVector &x) const override;
    
    bool hessian(const utopia::PetscVector &x, utopia::PetscMatrix &H) const override;

    bool value(const utopia::PetscVector &x, PetscScalar &value) const override;

    bool gradient(const utopia::PetscVector &x, utopia::PetscVector &g) const override;

    std::shared_ptr<dolfinx::fem::Function<T>> u();

private:
    std::unique_ptr<Impl> impl_;
};

#endif //UTOPIA_DOLFINX_FUNCTION_HPP
