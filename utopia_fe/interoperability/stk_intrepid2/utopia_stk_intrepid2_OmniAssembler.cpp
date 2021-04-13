#include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_ForcingFunction.hpp"
#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_LinearElasticity.hpp"
#include "utopia_intrepid2_VectorLaplaceOperator.hpp"
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_intrepid2.hpp"

#include <functional>

namespace utopia {
    namespace stk {

        // class AssemblerRegistry {
        // public:
        //     // using FEAssembler_t = utopia::FEAssembler<utopia::stk::FunctionSpace>;
        //     // using FEAssemblerPtr_t = std::shared_ptr<FEAssembler_t>;

        //     FEAssemblerPtr_t find_assembler(const std::string &name) const {
        //         auto it = assemblers_.find(name);
        //         if (it == assemblers_.end()) {
        //             // utopia::err() << "Unable to find assembler with name " << name << ".\n";
        //             // assert(false);
        //             return nullptr;
        //         } else {
        //             return it->second();
        //         }
        //     }

        //     template <typename Assembler_t>
        //     void register_assembler(const std::string &name) {
        //         assemblers_[name] = []() -> FEAssemblerPtr_t { return std::make_shared<Assembler_t>(); };
        //     }

        //     AssemblerRegistry() { register_assemblers(); }

        // private:
        //     std::map<std::string, std::function<FEAssemblerPtr_t()>> assemblers_;

        //     void register_assemblers() {
        //         register_assembler<Transport>("Transport");
        //         register_assembler<Mass>("Mass");
        //     }
        // };

        class OmniAssembler::Impl {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;

            template <class MaterialDescription>
            void init_material_assembler(const MaterialDescription &desc) {
                assemble_jacobian = [this, desc](const Vector &x, Matrix &mat) -> bool {
                    utopia::intrepid2::Assemble<MaterialDescription> assembler(desc, fe);
                    assembler.init();
                    local_to_global(*space, assembler.element_matrices(), ADD_MODE, mat);
                    return true;
                };

                assemble_material = [this, desc](const Vector &x, Matrix &mat, Vector &rhs) -> bool {
                    if (!assemble_jacobian(x, mat)) return false;
                    rhs = mat * x;
                    return true;
                };
            }

            template <class ForcingFunctionDescription>
            void init_forcing_function_assembler(const ForcingFunctionDescription &desc) {
                auto prev_fun = assemble_forcing_function;
                assemble_forcing_function = [this, prev_fun, desc](const Vector &x, Vector &rhs) -> bool {
                    if (prev_fun) {
                        // append
                        prev_fun(x, rhs);
                    }

                    utopia::intrepid2::Assemble<ForcingFunctionDescription> assembler(desc, fe);
                    assembler.init();
                    local_to_global(*space, assembler.element_vectors(), SUBTRACT_MODE, rhs);
                    return true;
                };
            }

            template <class ForcingFunctionDescription>
            void init_boundary_forcing_function_assembler(const std::string &boundary_name,
                                                          const ForcingFunctionDescription &desc) {
                auto prev_fun = assemble_forcing_function;
                assemble_forcing_function = [this, prev_fun, desc, boundary_name](const Vector &x,
                                                                                  Vector &rhs) -> bool {
                    if (prev_fun) {
                        // append
                        prev_fun(x, rhs);
                    }

                    FE boundary_fe;
                    create_fe_on_boundary(*space, boundary_fe, boundary_name, 2);

                    utopia::intrepid2::Assemble<ForcingFunctionDescription> assembler(desc,
                                                                                      utopia::make_ref(boundary_fe));
                    assembler.init();
                    side_local_to_global(*space, assembler.element_vectors(), SUBTRACT_MODE, rhs, boundary_name);
                    return true;
                };
            }

            std::function<bool(const Vector &x, Matrix &)> assemble_jacobian;
            std::function<bool(const Vector &x, Matrix &, Vector &)> assemble_material;
            std::function<bool(const Vector &x, Vector &)> assemble_forcing_function;

            std::shared_ptr<stk::FunctionSpace> space;
            std::shared_ptr<FE> fe;
            std::shared_ptr<Environment<stk::FunctionSpace>> env;
        };

        void OmniAssembler::set_environment(const std::shared_ptr<Environment<stk::FunctionSpace>> &env) {
            impl_->env = env;
        }

        OmniAssembler::OmniAssembler(const std::shared_ptr<stk::FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {
            impl_->space = space;
        }

        OmniAssembler::~OmniAssembler() = default;

        bool OmniAssembler::assemble(const Vector &x, Matrix &jacobian, Vector &fun) {
            if (!impl_->fe) {
                return false;
            }

            if (!impl_->assemble_material(x, jacobian, fun)) {
                return false;
            }

            if (impl_->assemble_forcing_function) {
                return impl_->assemble_forcing_function(x, fun);
            }

            return true;
        }

        bool OmniAssembler::assemble(const Vector &x, Matrix &jacobian) {
            if (!impl_->fe) {
                return false;
            }

            if (!impl_->assemble_jacobian(x, jacobian)) {
                return false;
            }

            return true;
        }

        bool OmniAssembler::assemble(const Vector &x, Vector &fun) {
            if (!impl_->fe) {
                return false;
            }

            assert(false && "IMPLEMENT ME");
            return false;
        }

        void OmniAssembler::read(Input &in) {
            // FIXME order must be guessed by discretization and material
            int quadrature_order = 2;
            in.get("quadrature_order", quadrature_order);
            impl_->fe = std::make_shared<Impl::FE>();
            create_fe(*impl_->space, *impl_->fe, quadrature_order);

            std::string material_type = "";

            in.get("material", [&](Input &in) { in.get("type", material_type); });

            const int spatial_dimension = impl_->space->mesh().spatial_dimension();

            // FIXME create a registry and decentralize material registration
            if (material_type == "LaplaceOperator") {
                LaplaceOperator<Scalar> material(1.0);
                in.get("material", material);
                impl_->init_material_assembler(material);

            } else if (material_type == "VectorLaplaceOperator") {
                switch (spatial_dimension) {
                    case 1: {
                        VectorLaplaceOperator<1, Scalar> material(1.0);
                        in.get("material", material);
                        impl_->init_material_assembler(material);
                        break;
                    }

                    case 2: {
                        VectorLaplaceOperator<2, Scalar> material(1.0);
                        in.get("material", material);
                        impl_->init_material_assembler(material);
                        break;
                    }

                    case 3: {
                        VectorLaplaceOperator<3, Scalar> material(1.0);
                        in.get("material", material);
                        impl_->init_material_assembler(material);
                        break;
                    }

                    default: {
                        assert(false);
                        break;
                    }
                }

            } else if (material_type == "LinearElasticity") {
                switch (spatial_dimension) {
                    case 1: {
                        LinearElasticity<1, Scalar> material(1.0, 1.0);
                        in.get("material", material);
                        impl_->init_material_assembler(material);
                        break;
                    }

                    case 2: {
                        LinearElasticity<2, Scalar> material(1.0, 1.0);
                        in.get("material", material);
                        impl_->init_material_assembler(material);
                        break;
                    }

                    case 3: {
                        LinearElasticity<3, Scalar> material(1.0, 1.0);
                        in.get("material", material);
                        impl_->init_material_assembler(material);
                        break;
                    }

                    default: {
                        assert(false);
                        break;
                    }
                }
            } else {
                if (material_type.empty()) {
                    utopia::err() << "[Error] Undefined material\n";
                } else {
                    utopia::err() << "[Error] Unsupported material " << material_type << '\n';
                }

                return;
            }

            in.get("forcing_functions", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    std::string forcing_function_type = "value";
                    std::string where = "domain";
                    std::string name;
                    int id = -1;

                    node.get("type", forcing_function_type);
                    node.get("where", where);
                    node.get("name", name);
                    node.get("id", id);

                    if (where == "surface") {
                        if (name.empty()) {
                            assert(id != -1);
                            name = "surface_" + std::to_string(id);
                        }

                        if (forcing_function_type == "value") {
                            ForcingFunction<Scalar> ff;
                            ff.read(node);
                            ff.n_components = impl_->space->n_var();
                            impl_->init_boundary_forcing_function_assembler(name, ff);
                        }

                    } else {
                        if (forcing_function_type == "value") {
                            ForcingFunction<Scalar> ff;
                            ff.read(node);
                            ff.n_components = impl_->space->n_var();
                            impl_->init_forcing_function_assembler(ff);
                        }
                    }
                });
            });
        }

    }  // namespace stk
}  // namespace utopia
