#include "utopia_intrepid2_OmniAssembler.hpp"

// utopia includes
#include "utopia_CreateFE.hpp"
#include "utopia_LocalToGlobal.hpp"

// utopia_intrepid2 includes
#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_ForcingFunction.hpp"
#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_LinearElasticity.hpp"
#include "utopia_intrepid2_Mass.hpp"
#include "utopia_intrepid2_NeoHookean.hpp"
#include "utopia_intrepid2_Transport.hpp"
#include "utopia_intrepid2_VectorLaplaceOperator.hpp"

#include <functional>

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace>
        class AssemblerRegistry {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using FE_t = utopia::intrepid2::FE<Scalar_t>;
            using FEAssembler_t = utopia::intrepid2::FEAssembler<Scalar_t>;
            using FEAssemblerPtr_t = std::shared_ptr<FEAssembler_t>;

            FEAssemblerPtr_t make_assembler(const std::shared_ptr<FE_t> &fe, Input &in) const {
                std::string type;
                in.get("type", type);

                if (has_variants(type)) {
                    type = create_variant_name(type, fe->spatial_dimension());
                }

                auto it = assemblers_.find(type);
                if (it == assemblers_.end()) {
                    assert(false);
                    return nullptr;
                } else {
                    return it->second(fe, in);
                }
            }

            template <typename MaterialDescription>
            void register_assembler(const std::string &type) {
                assemblers_[type] = [](const std::shared_ptr<FE_t> &fe, Input &in) -> FEAssemblerPtr_t {
                    MaterialDescription mat_desc;
                    mat_desc.read(in);
                    return std::make_shared<intrepid2::Assemble<MaterialDescription>>(fe, mat_desc);
                };
            }

            template <typename MaterialDescription>
            void register_assembler_variant(const std::string &type, const int spatial_dimension) {
                assemblers_[create_variant_name(type, spatial_dimension)] = [](const std::shared_ptr<FE_t> &fe,
                                                                               Input &in) -> FEAssemblerPtr_t {
                    MaterialDescription mat_desc;
                    mat_desc.read(in);
                    return std::make_shared<intrepid2::Assemble<MaterialDescription>>(fe, mat_desc);
                };

                has_ndim_variants.insert(type);
            }

            AssemblerRegistry() { register_assemblers(); }

            inline bool has_variants(const std::string &type) const {
                return has_ndim_variants.find(type) != has_ndim_variants.end();
            }

            static inline std::string create_variant_name(const std::string &type, int spatial_dimension) {
                return type + std::to_string(spatial_dimension);
            }

        private:
            std::map<std::string, std::function<FEAssemblerPtr_t(const std::shared_ptr<FE_t> &, Input &)>> assemblers_;
            std::set<std::string> has_ndim_variants;

            void register_assemblers() {
                register_assembler<utopia::Mass<Scalar_t>>("Mass");
                register_assembler<utopia::LaplaceOperator<Scalar_t>>("LaplaceOperator");
                register_assembler<utopia::ForcingFunction<Scalar_t>>("ForcingFunction");
                register_assembler<utopia::ForcingFunction<Scalar_t>>("ForcingFunction");

                register_assembler_variant<utopia::VectorLaplaceOperator<1, Scalar_t>>("VectorLaplaceOperator", 1);
                register_assembler_variant<utopia::VectorLaplaceOperator<2, Scalar_t>>("VectorLaplaceOperator", 2);
                register_assembler_variant<utopia::VectorLaplaceOperator<3, Scalar_t>>("VectorLaplaceOperator", 3);

                register_assembler_variant<utopia::LinearElasticity<1, Scalar_t>>("LinearElasticity", 1);
                register_assembler_variant<utopia::LinearElasticity<2, Scalar_t>>("LinearElasticity", 2);
                register_assembler_variant<utopia::LinearElasticity<3, Scalar_t>>("LinearElasticity", 3);

                register_assembler_variant<utopia::NeoHookean<1, Scalar_t>>("NeoHookean", 1);
                register_assembler_variant<utopia::NeoHookean<2, Scalar_t>>("NeoHookean", 2);
                register_assembler_variant<utopia::NeoHookean<3, Scalar_t>>("NeoHookean", 3);

                register_assembler_variant<utopia::Transport<1, typename FE_t::DynRankView>>("Transport", 1);
                register_assembler_variant<utopia::Transport<2, typename FE_t::DynRankView>>("Transport", 2);
                register_assembler_variant<utopia::Transport<3, typename FE_t::DynRankView>>("Transport", 3);
            }
        };

        template <class FunctionSpace>
        class OmniAssembler<FunctionSpace>::Impl {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using AssemblerRegistry = utopia::intrepid2::AssemblerRegistry<FunctionSpace>;

            // template <class MaterialDescription>
            // void init_material_assembler(const MaterialDescription &desc) {
            //     assemble_jacobian = [this, desc](const Vector &x, Matrix &mat) -> bool {
            //         utopia::intrepid2::Assemble<MaterialDescription> assembler(desc, fe);
            //         assembler.assemble();
            //         local_to_global(*space, assembler.data(), ADD_MODE, mat);
            //         return true;
            //     };

            //     assemble_material = [this, desc](const Vector &x, Matrix &mat, Vector &rhs) -> bool {
            //         if (!assemble_jacobian(x, mat)) return false;
            //         rhs = mat * x;
            //         return true;
            //     };
            // }

            template <class ForcingFunctionDescription>
            void init_forcing_function_assembler(const ForcingFunctionDescription &desc) {
                auto prev_fun = assemble_forcing_function;
                assemble_forcing_function = [this, prev_fun, desc](const Vector &x, Vector &rhs) -> bool {
                    if (prev_fun) {
                        // append
                        prev_fun(x, rhs);
                    }

                    utopia::intrepid2::Assemble<ForcingFunctionDescription> assembler(fe, desc);
                    assembler.assemble();
                    local_to_global(*space, assembler.vector_data(), SUBTRACT_MODE, rhs);
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

                    utopia::intrepid2::Assemble<ForcingFunctionDescription> assembler(utopia::make_ref(boundary_fe),
                                                                                      desc);
                    assembler.assemble();
                    side_local_to_global(*space, assembler.vector_data(), SUBTRACT_MODE, rhs, boundary_name);
                    return true;
                };
            }

            bool assemble_material(const Vector &x, Matrix &mat, Vector &vec) {
                if (!assemble_jacobian(x, mat)) return false;

                vec = mat * x;
                return true;
            }

            bool assemble_jacobian(const Vector &x, Matrix &mat) {
                if (assemblers.empty()) return false;
                assemblers[0]->ensure_matrix_accumulator();
                auto acc = assemblers[0]->matrix_accumulator();
                acc->zero();

                for (auto ass : assemblers) {
                    ass->set_matrix_accumulator(acc);
                    if (!ass->assemble()) {
                        assert(false);
                        return false;
                    }
                }

                local_to_global(*space, acc->data(), ADD_MODE, mat);
                return true;
            }

            std::function<bool(const Vector &x, Vector &)> assemble_forcing_function;

            std::vector<typename AssemblerRegistry::FEAssemblerPtr_t> assemblers;

            std::shared_ptr<FunctionSpace> space;
            std::shared_ptr<FE> fe;
            std::shared_ptr<Environment<FunctionSpace>> env;

            AssemblerRegistry registry;
        };

        template <class FunctionSpace>
        void OmniAssembler<FunctionSpace>::set_environment(const std::shared_ptr<Environment<FunctionSpace>> &env) {
            impl_->env = env;
        }

        template <class FunctionSpace>
        OmniAssembler<FunctionSpace>::OmniAssembler(const std::shared_ptr<FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {
            impl_->space = space;
        }

        template <class FunctionSpace>
        OmniAssembler<FunctionSpace>::~OmniAssembler() = default;

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::assemble(const Vector &x, Matrix &jacobian, Vector &fun) {
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

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::assemble(const Vector &x, Matrix &jacobian) {
            if (!impl_->fe) {
                return false;
            }

            if (!impl_->assemble_jacobian(x, jacobian)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::assemble(const Vector &x, Vector &fun) {
            if (!impl_->fe) {
                return false;
            }

            assert(false && "IMPLEMENT ME");
            return false;
        }

        template <class FunctionSpace>
        void OmniAssembler<FunctionSpace>::read(Input &in) {
            // FIXME order must be guessed by discretization and material
            int quadrature_order = 2;
            in.get("quadrature_order", quadrature_order);
            impl_->fe = std::make_shared<typename Impl::FE>();
            create_fe(*impl_->space, *impl_->fe, quadrature_order);

            in.get("material", [this](Input &node) {
                auto assembler = impl_->registry.make_assembler(impl_->fe, node);

                if (assembler) {
                    impl_->assemblers.push_back(assembler);
                } else {
                    assert(false && "Should not come here");
                    Utopia::Abort();
                }
            });

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

    }  // namespace intrepid2
}  // namespace utopia