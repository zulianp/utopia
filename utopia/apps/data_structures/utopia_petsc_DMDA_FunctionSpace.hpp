#ifndef UTOPIA_PETSC_DMDA_FUNCTION_SPACE_HPP
#define UTOPIA_PETSC_DMDA_FUNCTION_SPACE_HPP

#include "utopia_FunctionSpace.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_MultiVariateElement.hpp"

//TOBEREMOVED
#include "utopia_petsc_dma_FunctionSpace.hpp"


namespace utopia {

    template<class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    class FunctionSpace< PetscDMDA<Point_, IntArray_>, NComponents_, UniVarElem_> : public Configurable {
    public:
        //concrete types
        using Vector = utopia::PetscVector;
        using Matrix = utopia::PetscMatrix;
        using Comm   = utopia::PetscCommunicator;

        //from template arg list
        using Point = Point_;
        using Mesh  = utopia::PetscDMDA<Point_, IntArray_>;
        using Shape = UniVarElem_;
        using Elem  = MultiVariateElem<UniVarElem_, NComponents_>;

        using NodeIndex = typename Mesh::NodeIndex;

        using MemType  = typename Elem::MemType;
        using Scalar   = typename Mesh::Scalar;
        using SizeType = typename Mesh::SizeType;


        using ViewDevice = FunctionSpace;
        using Device = typename Mesh::Device;


        using DirichletBC = utopia::DirichletBoundaryCondition<FunctionSpace>;
        using DofMapping  = utopia::DofMapping<Mesh, UniVarElem_, NComponents_>;
        static const int NDofs = DofMapping::NDofs;

    };

    template<class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    using PetscDMDAFunctionSpace = FunctionSpace< PetscDMDA<Point_, IntArray_>, NComponents_, UniVarElem_>;

}

#endif //UTOPIA_PETSC_DMDA_FUNCTION_SPACE_HPP
