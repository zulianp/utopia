#include "utopia_libmesh_AssemblyContext.hpp"

namespace utopia {

    void LibMeshAssemblyContext::init_tensor(Vector &v, const bool reset)
    {
        const auto nt = n_shape_functions();
        auto s = size(v);

        if(reset || nt != s.get(0)) {
            v = zeros(nt);
        }
    }

    void LibMeshAssemblyContext::init_tensor(Matrix &v, const bool reset)
    {
        auto s = size(v);
        const auto n_test = n_shape_functions();
        const auto n_trial = n_test;

        if(reset || n_test != s.get(0) || n_trial != s.get(1)) {
            v = zeros(n_test, n_trial);
        }
    }
}
