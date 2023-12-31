#include "utopia_FractureFlow.hpp"

namespace utopia {

    FractureFlow::FractureFlow(libMesh::Parallel::Communicator &comm) : mesh(comm), space(make_ref(mesh)) {}

    void FractureFlow::read(Input &is) {
        try {
            is.get("mesh", mesh);
            is.get("space", space);

            auto grid_sampler = std::make_shared<UIScalarSampler<double>>();
            is.get("sampler", *grid_sampler);

            if (!grid_sampler->empty()) {
                sampler = grid_sampler;
            } else {
                auto subdomain_fun = utopia::make_unique<UISubdomainFunction<double>>();

                is.get("diffusivity-blocks", *subdomain_fun);

                if (!subdomain_fun->good()) {
                    sampler = std::make_shared<UIConstantFunction<double>>(1.);
                } else {
                    if (!subdomain_fun->has_default()) {
                        subdomain_fun->set_default(utopia::make_unique<UIConstantFunction<double>>(1.));
                    }

                    sampler = std::move(subdomain_fun);
                }
            }

            forcing_function = std::make_shared<UIForcingFunction<FunctionSpaceT, UVector>>(space.subspace(0));
            is.get("forcing-function", *forcing_function);

            is.get("weak-bc", [this](Input &in) {
                weak_BC_ = std::make_shared<WeakDirichletBoundaryConditions<FunctionSpaceT, USparseMatrix, UVector>>(
                    space.subspace(0));
                weak_BC_->read(in);
            });

            // material parameters
            double diffusivity = 1.;
            double diffusivities[3] = {1., 1., 1.};

            is.get("diffusivity", diffusivity);
            is.get("diffusivity-x", diffusivities[0]);
            is.get("diffusivity-y", diffusivities[1]);
            is.get("diffusivity-z", diffusivities[2]);

            int dim = space.subspace(0).mesh().spatial_dimension();

            diffusion_tensor = identity(dim, dim);

            {
                Write<ElementMatrix> w(diffusion_tensor);
                for (int i = 0; i < dim; ++i) {
                    diffusion_tensor.set(i, i, diffusivities[i] * diffusivity);
                }
            }

            is.get("post-processors", [this](Input &in) {
                in.get_all([this](Input &in) {
                    std::string type;
                    in.get("type", type);

                    if (type == "flux") {
                        auto flux = std::make_shared<FluxPostProcessor<FunctionSpaceT, UVector>>();

                        flux->sampler(sampler);
                        flux->diffusion_tensor(diffusion_tensor);

                        flux->read(in);

                        post_processors_.push_back(flux);
                    } else if (type == "avg") {
                        auto flux = std::make_shared<AverageHeadPostProcessor<FunctionSpaceT, UVector>>();

                        flux->read(in);

                        post_processors_.push_back(flux);
                    }
                });
            });

        } catch (const std::exception &ex) {
            std::cerr << ex.what() << std::endl;
            assert(false);
        }
    }

    void FractureFlow::post_process(const UVector &sol) {
        for (auto pp : post_processors_) {
            pp->apply(space.space()[0], sol);
            pp->describe();
        }
    }

    void FractureFlow::post_process(LibMeshFunctionSpace &space,
                                    const UVector &pressure,
                                    const UVector &concentration) {
        for (auto pp : post_processors_) {
            pp->apply(space, pressure, concentration);
            pp->describe();
        }
    }

    void FractureFlow::export_post_process() {
        for (auto pp : post_processors_) {
            pp->export_values();
        }
    }

    void FractureFlow::apply_weak_BC(USparseMatrix &A, UVector &b) const {
        if (weak_BC_) {
            weak_BC_->apply(A, b);
        }
    }

    bool FractureFlow::empty() const { return mesh.empty(); }

    void FractureFlow::describe(std::ostream &os) const {
        os << "-----------------------------------\n";
        // mesh.describe(os);
        // space.describe(os);
        // forcing_function.describe(os);
        os << "permeability: " << std::endl;

        {
            int dim = size(diffusion_tensor).get(0);
            Read<ElementMatrix> w(diffusion_tensor);
            for (int i = 0; i < dim; ++i) {
                os << diffusion_tensor.get(i, i) << " ";
            }
        }

        os << "\n";
        os << "-----------------------------------\n";
    }

    FractureFlowAuxSystem::FractureFlowAuxSystem(LibMeshFunctionSpace &V, const std::string &name)
        : aux_(V.equation_systems().add_system<libMesh::LinearImplicitSystem>("aux")) {
        var_nums_.push_back(aux_.add_variable(name, libMesh::Order(V.order(0)), libMesh::LAGRANGE));
        aux_.init();
    }

    void FractureFlowAuxSystem::sample(const std::shared_ptr<UIFunction<double>> &sampler) {
        LibMeshFunctionSpace V_aperture(aux_, var_nums_[0]);
        V_aperture.initialize();

        auto &dof_map = V_aperture.dof_map();

        UIndexArray ghost_nodes;
        convert(dof_map.get_send_list(), ghost_nodes);
        sampled = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
        // sampled = local_zeros(dof_map.n_local_dofs());

        auto constant_sampler = std::dynamic_pointer_cast<UIConstantFunction<double>>(sampler);

        if (constant_sampler) {
            sampled.set(constant_sampler->value());
        } else {
            auto u = trial(V_aperture);
            auto v = test(V_aperture);

            auto lform = inner(ctx_fun(sampler), v) * dX;

            UIndexArray ghost_nodes;
            convert(dof_map.get_send_list(), ghost_nodes);
            UVector aperture_h = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
            utopia::assemble(lform, aperture_h);

            USparseMatrix mass_mat;
            utopia::assemble(inner(u, v) * dX, mass_mat);
            UVector d_inv = 1. / sum(mass_mat, 1);
            sampled = e_mul(d_inv, aperture_h);
        }

        utopia::convert(sampled, *aux_.solution);
        aux_.solution->close();
    }

}  // namespace utopia
