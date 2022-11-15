#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

// Stk includes
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

// utopia/kokkos includes
#include "utopia_kokkos_L2Norm.hpp"
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_LinearElasticity.hpp"
#include "utopia_kokkos_NeoHookean.hpp"
#include "utopia_kokkos_VectorLaplaceOperator.hpp"

// utopia/intrepid2 includes
#include "utopia_intrepid2_ShellTools.hpp"

// Interop includes
#include "utopia_stk_intrepid2.hpp"

// FIXME: This is the last include because the operator files are not yet in the correct place
#include "utopia_SpaceAndFETest.hpp"

#include "utopia_ScalarProductTest.hpp"

#include "utopia_kokkos_LaplaceOperator_new.hpp"

#include "utopia_stk_intrepid2_Discretization.hpp"

using namespace utopia;

using StkScalar = Traits<utopia::stk::FunctionSpace>::Scalar;
void interop_stk_intrepid2() { SpaceAndFETest<utopia::stk::FunctionSpace, utopia::intrepid2::FE<StkScalar>>().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

void scalar_product_stk_intrepid2() { ScalarProductTest<utopia::stk::FunctionSpace>().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(scalar_product_stk_intrepid2);

void new_assembler_test() {
    using FS_t = utopia::stk::FunctionSpace;
    using FE_t = utopia::intrepid2::FE<double>;
    using Matrix_t = Traits<FS_t>::Matrix;
    using Vector_t = Traits<FS_t>::Vector;
    using Assembler_t = utopia::kokkos::FEAssembler<FS_t, FE_t>;
    using Discretization_t = utopia::Discretization<FS_t, FE_t>;

    UnitCubeSpaceAndFETest<FS_t, FE_t> test;
    auto params = test.cube_space_param(1);

    FS_t space;
    space.read(params);

    test.add_cube_bc(space, 1);

    auto fe_ptr = std::make_shared<FE_t>();
    create_fe(space, *fe_ptr, 2);

    auto discretization = std::make_shared<Discretization_t>(make_ref(space), fe_ptr);
    // FIXME remove args once we remove the superset features
    auto assembler = std::make_shared<Assembler_t>(fe_ptr);
    assembler->set_discretization(discretization);

    utopia::kokkos::LaplaceOperatorNew<FS_t, FE_t, Assembler_t> lapl;
    lapl.set_assembler(assembler);

    Matrix_t mat;
    space.create_matrix(mat);

    Vector_t x;
    space.create_vector(x);
    x.set(0.0);

    lapl.hessian(x, mat);

    // disp(mat);
}

UTOPIA_REGISTER_TEST_FUNCTION(new_assembler_test);

#endif  // UTOPIA_WITH_INTREPID2
