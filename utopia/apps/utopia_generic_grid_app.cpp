
#include "utopia_Base.hpp"

#include "utopia_petsc_Base.hpp"

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)

#include "petscfe.h"
#include "utopia_Algorithms.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_Rename.hpp"
#include "utopia_StructuredGrid.hpp"
#include "utopia_Views.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DMPlex.hpp"
#include "utopia_petsc_DMPlex_FunctionSpace.hpp"
#include "utopia_ui.hpp"

namespace utopia {
    using V = utopia::StaticVector<PetscScalar, 2>;
    using I = utopia::ArrayView<PetscInt, 2>;
    using VA = utopia::ArrayView<PetscScalar, 2>;
    template class StructuredGrid<V, I>;
    template class PetscDMDA<V, I>;
    template class PetscDMPlex<V, I>;

    // FIXME
    void generic_grid_test(Input &in) {
        // COMPILE-TIME BEGIN (everything is evaualted at compile time)
        constexpr StructuredGrid<V, I> g({10, 10}, {0, 0}, {10, 10}, {0, 0}, {10, 10}, VA({0.0, 0.0}), VA({1.0, 1.0}));

        static_assert(g.dim() == 2, "dim must be constexpr");
        static_assert(g.n_nodes() == 100, "n_nodes has to be computable at compile time");
        static_assert(g.n_elements() == 81, "wrong number of elements");

        constexpr bool is_b = g.is_node_on_boundary(0);
        static_assert(is_b, "node 0 must be on boundary");
        // constexpr bool is_b = g.is_node_on_boundary_local_no_ghost(0, SideSet::left());

        // for some reaons clang is able to do this at compile time but not gcc
        /*constexpr*/ PetscScalar meas = g.measure();
        // COMPILE-TIME END
        assert(device::approxeq(meas, 1.0, 1e-8));
    }

    UTOPIA_REGISTER_APP(generic_grid_test);

    void vec_grid_test(Input &in) {
        // run-time sizes real data
        std::vector<PetscInt> zeros_v = {0, 0};
        std::vector<PetscInt> dims_v = {10, 10};
        std::vector<PetscScalar> box_min_v = {0.0f, 0.0f};
        std::vector<PetscScalar> box_max_v = {1.0f, 1.0f};

        // views on data
        ArrayView<PetscInt> zeros(&zeros_v[0], zeros_v.size());
        ArrayView<PetscInt> dims(&dims_v[0], dims_v.size());
        ArrayView<PetscScalar> box_min(&box_min_v[0], box_min_v.size());
        ArrayView<PetscScalar> box_max(&box_max_v[0], box_max_v.size());

        StructuredGrid<VectorView<ArrayView<PetscScalar>>, ArrayView<PetscInt>> g(
            dims, zeros, dims, zeros, dims, box_min, box_max);

        assert(g.dim() == 2);
        assert(g.n_nodes() == 100);
        assert(g.n_elements() == 81);
        PetscScalar meas = g.measure();
        bool is_b = g.is_node_on_boundary(0);
        assert(is_b);

        assert(device::approxeq(meas, 1.0, 1e-8));
    }

    UTOPIA_REGISTER_APP(vec_grid_test);

    void dmda_test(Input &in) {
        PetscCommunicator comm;
        PetscDMDA<V, I> dmda(comm);

        dmda.read(in);

        PetscVector v;
        dmda.create_vector(v);
        v.set(1.0);

        dmda.write("prova.vtr", v);
    }

    UTOPIA_REGISTER_APP(dmda_test);

    void dmplex_test(Input &in) {
        PetscCommunicator comm;
        MPITimeStatistics stats(comm);

        stats.start();

        PetscDMPlex<V, I> dmplex(comm);
        dmplex.read(in);

        stats.stop_collect_and_restart("dmplex-setup");

        dmplex.write("mesh.vtu");

        stats.stop_and_collect("vtu-write");
        stats.describe(std::cout);

        PetscInt num_fields = 1;
        PetscInt dim = dmplex.dim(), qorder = PETSC_DEFAULT;
        in.get("dim", dim);
        in.get("qorder", qorder);

        DMSetNumFields(dmplex.raw_type(), num_fields);

        // Important to specifiy for all num_fields
        PetscInt num_comp[1] = {1};

        // Important to specifiy all mesh dim 0, 1, 2, 3
        PetscInt num_dofs[4] = {1, 1, 0, 0};

        PetscSection section = nullptr;
        DMPlexCreateSection(
            dmplex.raw_type(), nullptr, num_comp, num_dofs, 0, nullptr, nullptr, nullptr, nullptr, &section);

        PetscSectionSetFieldName(section, 0, "u");
        DMSetLocalSection(dmplex.raw_type(), section);
        PetscSectionDestroy(&section);
        DMSetUp(dmplex.raw_type());

        PetscVector v;
        dmplex.create_vector(v);
        v.set(1.0);

        v.comm().root_print("dofs = " + std::to_string(v.size()));

        utopia::rename("X", v);
        dmplex.write("prova.vtu", v);

        PetscReal ref_coords[3] = {0.0, 0.0, 0.0};
        PetscReal real_coords[3] = {0.0, 0.0, 0.0};
        DMPlexReferenceToCoordinates(dmplex.raw_type(), 0, 1, ref_coords, real_coords);

        // ref element is center in 0 with -1, 1 range

        std::cout << real_coords[0] << " " << real_coords[1] << " " << real_coords[2] << std::endl;

        PetscReal vert[3], J[3 * 3], invJ[3 * 3], detJ;
        DMPlexComputeCellGeometryFEM(dmplex.raw_type(), 0, nullptr, vert, J, invJ, &detJ);

        std::cout << "detJ: " << detJ << std::endl;

        DMPlexComputeCellGeometryAffineFEM(dmplex.raw_type(), 0, vert, J, invJ, &detJ);

        std::cout << "detJ: " << detJ << std::endl;

        comm.barrier();

        // PetscVector coords;
        // coords.destroy();
        // DMGetCoordinates(dmplex.raw_type(), &coords.raw_type());

        // DMGetCoordinatesLocal(dmplex.raw_type(), &coords.raw_type());

        // if (comm.rank() == 0) disp(coords);
        // comm.barrier();
        // if (comm.rank() == 1) disp(coords);

        // coords.raw_type() = nullptr;

        // VecView()

        PetscInt num_points = 0;
        PetscInt *points = nullptr;

        PetscInt cell_num = 0;
        in.get("cell_num", cell_num);
        DMPlexGetTransitiveClosure(dmplex.raw_type(), cell_num, PETSC_TRUE, &num_points, &points);

        std::cout << "num_points: " << num_points << std::endl;

        for (PetscInt i = 1; i < num_points; ++i) {
            std::cout << points[i * 2] << " ";
        }

        std::cout << std::endl;

        DMPlexRestoreTransitiveClosure(dmplex.raw_type(), 0, PETSC_TRUE, &num_points, &points);

        PetscInt s, e;

        PetscInt point = 4;
        in.get("point", point);
        DMPlexGetPointLocal(dmplex.raw_type(), point, &s, &e);

        std::cout << "[" << s << ", " << e << ")" << std::endl;
    }

    UTOPIA_REGISTER_APP(dmplex_test);
}  // namespace utopia

#endif  // UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)
