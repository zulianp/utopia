
#include "utopia_Base.hpp"


#include "utopia_StructuredGrid.hpp"
#include "utopia_ui.hpp"
#include "utopia_Views.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMPlex.hpp"
#include "utopia_Rename.hpp"

namespace utopia {
    using V = utopia::StaticVector<PetscScalar, 2>;
    using I = utopia::ArrayView<PetscInt, 2>;
    using VA = utopia::ArrayView<PetscScalar, 2>;
    template class StructuredGrid<V, I>;
    template class PetscDMDA<V, I>;
    template class PetscDMPlex<V, I>;

    //FIXME
    void generic_grid_test(Input &in)
    {
        //COMPILE-TIME BEGIN (everything is evaualted at compile time)
        constexpr StructuredGrid<V, I> g(
            {10, 10},
            {0, 0}, {10, 10},
            {0, 0}, {10, 10},
            VA({0.0, 0.0}),
            VA({1.0, 1.0})
        );

        static_assert(g.dim() == 2, "dim must be constexpr");
        static_assert(g.n_nodes() == 100, "n_nodes has to be computable at compile time");
        static_assert(g.n_elements() == 81, "wrong number of elements");

        constexpr bool is_b = g.is_node_on_boundary(0);
        static_assert(is_b, "node 0 must be on boundary");
        // constexpr bool is_b = g.is_node_on_boundary_local_no_ghost(0, SideSet::left());

        //for some reaons clang is able to do this at compile time but not gcc
        /*constexpr*/ PetscScalar meas = g.measure();
        //COMPILE-TIME END
        assert(device::approxeq(meas, 1.0, 1e-8));

    }

    UTOPIA_REGISTER_APP(generic_grid_test);


    void vec_grid_test(Input &in)
    {
        //run-time sizes real data
        std::vector<PetscInt>   zeros_v  = {0, 0};
        std::vector<PetscInt>   dims_v   = {10, 10};
        std::vector<PetscScalar> box_min_v = {0.0f, 0.0f};
        std::vector<PetscScalar> box_max_v = {1.0f, 1.0f};


        //views on data
        ArrayView<PetscInt>   zeros(&zeros_v[0], zeros_v.size());
        ArrayView<PetscInt>   dims(&dims_v[0], dims_v.size());
        ArrayView<PetscScalar> box_min(&box_min_v[0], box_min_v.size());
        ArrayView<PetscScalar> box_max(&box_max_v[0], box_max_v.size());

        StructuredGrid<VectorView<ArrayView<PetscScalar>>, ArrayView<PetscInt>> g(
            dims,
            zeros, dims,
            zeros, dims,
            box_min,
            box_max
        );

        assert(g.dim() == 2);
        assert(g.n_nodes() == 100);
        assert(g.n_elements() == 81);
        PetscScalar meas = g.measure();
        bool is_b = g.is_node_on_boundary(0);
        assert(is_b);

        assert(device::approxeq(meas, 1.0, 1e-8));

    }

    UTOPIA_REGISTER_APP(vec_grid_test);


    void dmda_test(Input &in)
    {
        PetscCommunicator comm;
        PetscDMDA<V, I> dmda(comm);

        dmda.read(in);

        PetscVector v;
        dmda.create_vector(v);
        v.set(1.0);

        dmda.write("prova.vtr", v);
    }

    UTOPIA_REGISTER_APP(dmda_test);

    void dmplex_test(Input &in)
    {
        PetscCommunicator comm;
        PetscDMPlex<V, I> dmplex(comm);

        // dmplex.read(in);
        dmplex.read("../../utopia_fe/data/contact/contact_squares_tri.e");


        // DMLabel *label = nullptr;
        PetscInt num_fields = 2;
        // PetscInt num_comp[1] = {1};
        // PetscInt num_dof[1]  = {1};
        // PetscInt num_bc = 0;
        // PetscInt *bc_field = nullptr;
        // IS *bc_comps = nullptr,* bc_points = nullptr, perm = nullptr;

          // DMGetStratumIS(dm, "marker", 1, &bcPointIS[0]);

        DMSetNumFields(dmplex.raw_type(),  num_fields);

        PetscFE   fe[2];
        PetscFECreateDefault(dmplex.comm().get(), 2,1,PETSC_TRUE,nullptr,PETSC_DEFAULT,&fe[0]);
        PetscFESetName(fe[0], "u");
        DMSetField(dmplex.raw_type(),0,nullptr,(PetscObject)fe[0]);

        PetscFECreateDefault(dmplex.comm().get(), 2,1,PETSC_TRUE,nullptr,PETSC_DEFAULT,&fe[1]);
        PetscFESetName(fe[1], "v");

        DMSetField(dmplex.raw_type(),1,nullptr,(PetscObject)fe[1]);
        DMCreateDS(dmplex.raw_type());

        PetscFEDestroy(&fe[0]);
        PetscFEDestroy(&fe[1]);


        // DMPlexCreateClosureIndex(dmplex.raw_type(), nullptr);

        // DMSetUp(dmplex.raw_type());

        // PetscSection section;
        // DMPlexCreateSection(dmplex.raw_type(), label, num_comp, num_dof, num_bc, bc_field, bc_comps, bc_points, perm, &section);
        // PetscSectionSetFieldName(section, 0, "u");
        // // DMSetLocalSection(dmplex.raw_type(), section);
        // DMSetGlobalSection(dmplex.raw_type(), section);

        // dmplex.create_section()

        PetscVector v;
        dmplex.create_vector(v);
        v.set(1.0);

        utopia::rename("X", v);

        // dmplex.write("prova.vtu");//, v);
        dmplex.write("prova.vtu", v);

        // PetscSectionDestroy(&section);
    }

    UTOPIA_REGISTER_APP(dmplex_test);
}

