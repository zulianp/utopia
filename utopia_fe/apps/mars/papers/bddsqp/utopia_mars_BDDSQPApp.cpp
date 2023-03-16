#include "utopia_BDDSQPApp.hpp"

#include "utopia_Main.hpp"

#include "utopia_mars.hpp"

#include "utopia_mars_Discretization.hpp"
#include "utopia_mars_Material.hpp"

void mars_bddsqp(utopia::Input &in) {
    utopia::BDDSQPApp<utopia::mars::FunctionSpace, utopia::kokkos::UniformFE<double>> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(mars_bddsqp);
