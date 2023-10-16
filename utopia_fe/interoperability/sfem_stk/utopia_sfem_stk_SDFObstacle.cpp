#include "utopia_SDFObstacle_impl.hpp"

#include "utopia_ContactFactory.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include "utopia_moonolith_stk.hpp"

#include "utopia_sfem_stk_ExtractSurface.hpp"

namespace utopia {
    template class SDFObstacle<utopia::stk::FunctionSpace>;
}

using FunctionSpace_t = utopia::stk::FunctionSpace;
using SDFObstacle_t = utopia::SDFObstacle<utopia::stk::FunctionSpace>;

UTOPIA_CONTACT_REGISTER(FunctionSpace_t, SDFObstacle_t, "sdf");
