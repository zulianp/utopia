#include "utopia_plugin_Function_impl.hpp"

#include "utopia_Communicator.hpp"
#include "utopia_IOStream.hpp"
#include "utopia_Instance.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Path.hpp"

#include <dlfcn.h>

#define CHECK_REQUIRED(plugin_fun, name)                                     \
    do {                                                                     \
        if (!plugin_fun) {                                                   \
            if (utopia::mpi_world_rank() == 0) {                             \
                utopia::err() << "Could not find function " << name << '\n'; \
                utopia::Utopia::Abort();                                     \
            }                                                                \
        }                                                                    \
    } while (0)

static void check_error() {
    auto msg = dlerror();

    if (msg) {
        if (utopia::mpi_world_rank() == 0) {
            utopia::err() << "Error when reading plugin \n";
            utopia::err() << "msg: " << msg << '\n';
        }
    }
}

static void *load_function(void *handle, const char *name) {
    auto ret = dlsym(handle, name);
    check_error();
    CHECK_REQUIRED(ret, name);
    return ret;
}

namespace utopia {
    void PluginFunctionImpl::read(Input &in) {
        utopia::Path path;
        in.require("path", path);

        handle = dlopen(path.c_str(), RTLD_NOW);
        check_error();

        if (!handle) {
            if (mpi_world_rank() == 0) {
                utopia::err() << "Unable to load plugin " << path << '\n';
            }

            Utopia::Abort();
        }

        init_ = (utopia_plugin_Function_init_t)load_function(handle, "utopia_plugin_Function_init");

        create_crs_graph_ =
            (utopia_plugin_Function_create_crs_graph_t)load_function(handle, "utopia_plugin_Function_create_crs_graph");

        create_vector_ =
            (utopia_plugin_Function_create_vector_t)load_function(handle, "utopia_plugin_Function_create_vector");

        destroy_vector_ =
            (utopia_plugin_Function_destroy_vector_t)load_function(handle, "utopia_plugin_Function_destroy_vector");

        value_ = (utopia_plugin_Function_value_t)load_function(handle, "utopia_plugin_Function_value");

        gradient_ = (utopia_plugin_Function_gradient_t)load_function(handle, "utopia_plugin_Function_gradient");

        hessian_crs_ =

            (utopia_plugin_Function_hessian_crs_t)load_function(handle, "utopia_plugin_Function_hessian_crs");
        apply_ = (utopia_plugin_Function_apply_t)load_function(handle, "utopia_plugin_Function_apply");

        apply_constraints_ = (utopia_plugin_Function_apply_constraints_t)load_function(
            handle, "utopia_plugin_Function_apply_constraints");

        apply_zero_constraints_ = (utopia_plugin_Function_apply_zero_constraints_t)load_function(
            handle, "utopia_plugin_Function_apply_zero_constraints");

        copy_constrained_dofs_ = (utopia_plugin_Function_copy_constrained_dofs_t)load_function(
            handle, "utopia_plugin_Function_copy_constrained_dofs");

        destroy_ = (utopia_plugin_Function_destroy_t)load_function(handle, "utopia_plugin_Function_destroy");
    }

    PluginFunctionImpl::~PluginFunctionImpl() {
        if (this->destroy_) {
            this->destroy_(&info);
        }

        dlclose(handle);
    }

    void PluginFunctionImpl::initialize(const Communicator &comm) {
        info.comm = comm.raw_comm();
        int err = this->init_(&info);
        if (err != UTOPIA_PLUGIN_SUCCESS) {
            if (mpi_world_rank() == 0) {
                utopia::err() << "Plugin return error code " << err << '\n';
            }

            Utopia::Abort();
        }
    }
}  // namespace utopia
