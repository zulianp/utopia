
#ifndef UTOPIA_MARS_MESH_IO_HPP
#define UTOPIA_MARS_MESH_IO_HPP

// #include "utopia_Input.hpp"
// #include "utopia_Path.hpp"
// #include "utopia_Traits.hpp"

// #include "utopia_fe_Core.hpp"

// #include "utopia_mars_ForwardDeclarations.hpp"
// #include "utopia_mars_Mesh.hpp"

// #include <memory>
// #include <string>

#include "mars_config.hpp"

#ifdef MARS_ENABLE_ADIOS2
#include "mars_adios2_IO.hpp"
#else
#ifdef MARS_ENABLE_VTK
#include "mars_vtk_IO.hpp"
#endif
#endif

namespace utopia {
    namespace mars {

#ifdef MARS_ENABLE_ADIOS2
#define MARS_WITH_WITH_IO
        template <class... Args>
        using MarsIOImpl = ::mars::adios2::IO<Args...>;
#else
#ifdef MARS_ENABLE_VTK
#define MARS_WITH_WITH_IO
        template <class... Args>
        using MarsIOImpl = ::mars::vtk::IO<Args...>;
#endif
#endif

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_MESH_IO_HPP
