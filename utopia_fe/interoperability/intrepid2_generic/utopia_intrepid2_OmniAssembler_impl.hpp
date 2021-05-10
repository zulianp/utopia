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
#include "utopia_intrepid2_Residual.hpp"
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
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar_t>;
            using AssemblerRegistry = utopia::intrepid2::AssemblerRegistry<FunctionSpace>;
            using FEAssembler_t = utopia::intrepid2::FEAssembler<Scalar_t>;
            using FEAssemblerPtr_t = std::shared_ptr<FEAssembler_t>;
            using TensorAccumulator_t = typename FEAssembler_t::TensorAccumulator;

            class PartAssembler {
            public:
                virtual ~PartAssembler() {}

                std::string name;
                std::shared_ptr<FE> fe;
                std::vector<FEAssemblerPtr_t> assemblers;

                std::shared_ptr<TensorAccumulator_t> matrix_accumulator;
                std::shared_ptr<TensorAccumulator_t> vector_accumulator;
                std::shared_ptr<TensorAccumulator_t> scalar_accumulator;

                void ensure_accumulators() {
                    for (auto &a_ptr : assemblers) {
                        if (a_ptr->is_matrix() && !matrix_accumulator) {
                            a_ptr->ensure_matrix_accumulator();
                            matrix_accumulator = a_ptr->matrix_accumulator();
                        } else {
                            a_ptr->set_matrix_accumulator(matrix_accumulator);
                        }

                        if (a_ptr->is_vector() && !vector_accumulator) {
                            a_ptr->ensure_vector_accumulator();
                            vector_accumulator = a_ptr->vector_accumulator();
                        } else {
                            a_ptr->set_vector_accumulator(vector_accumulator);
                        }

                        if (a_ptr->is_scalar() && !scalar_accumulator) {
                            a_ptr->ensure_scalar_accumulator();
                            scalar_accumulator = a_ptr->scalar_accumulator();
                        } else {
                            a_ptr->set_scalar_accumulator(scalar_accumulator);
                        }
                    }
                }

                void ensure_matrix_accumulators() {
                    for (auto &a_ptr : assemblers) {
                        if (a_ptr->is_matrix() && !matrix_accumulator) {
                            a_ptr->ensure_matrix_accumulator();
                            matrix_accumulator = a_ptr->matrix_accumulator();
                        } else {
                            a_ptr->set_matrix_accumulator(matrix_accumulator);
                        }
                    }
                }

                void ensure_vector_accumulators() {
                    for (auto &a_ptr : assemblers) {
                        if (a_ptr->is_vector() && !vector_accumulator) {
                            a_ptr->ensure_vector_accumulator();
                            vector_accumulator = a_ptr->vector_accumulator();
                        } else {
                            a_ptr->set_vector_accumulator(vector_accumulator);
                        }
                    }
                }

                void ensure_scalar_accumulators() {
                    for (auto &a_ptr : assemblers) {
                        if (a_ptr->is_scalar() && !scalar_accumulator) {
                            a_ptr->ensure_scalar_accumulator();
                            scalar_accumulator = a_ptr->scalar_accumulator();
                        } else {
                            a_ptr->set_scalar_accumulator(scalar_accumulator);
                        }
                    }
                }

                void ensure_vector_accumulator() {
                    if (assemblers.empty()) {
                        assert(false);
                        return;
                    }

                    if (!vector_accumulator) {
                        (*assemblers.begin())->ensure_vector_accumulator();
                        vector_accumulator = (*assemblers.begin())->vector_accumulator();
                    }
                }

                void zero_scalar_accumulators() {
                    if (has_scalar()) {
                        scalar_accumulator->zero();
                    }
                }

                void zero_matrix_accumulators() {
                    if (has_matrix()) {
                        matrix_accumulator->zero();
                    }
                }

                void zero_vector_accumulators() {
                    if (has_vector()) {
                        vector_accumulator->zero();
                    }
                }

                void zero_accumulators() {
                    zero_scalar_accumulators();
                    zero_vector_accumulators();
                    zero_matrix_accumulators();
                }

                inline bool has_matrix() const { return static_cast<bool>(matrix_accumulator); }
                inline bool has_vector() const { return static_cast<bool>(vector_accumulator); }
                inline bool has_scalar() const { return static_cast<bool>(scalar_accumulator); }
            };

            class DomainAssembler : public PartAssembler {};

            class BoundaryAssembler : public PartAssembler {};

            void ensure_accumulators() {
                domain.ensure_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.ensure_accumulators();
                }
            }

            void ensure_matrix_accumulators() {
                domain.ensure_matrix_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.ensure_matrix_accumulators();
                }
            }

            void ensure_vector_accumulators() {
                domain.ensure_vector_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.ensure_vector_accumulators();
                }
            }

            void ensure_scalar_accumulators() {
                domain.ensure_scalar_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.ensure_scalar_accumulators();
                }
            }

            void zero_accumulators() {
                domain.zero_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.zero_accumulators();
                }
            }

            void zero_scalar_accumulators() {
                domain.zero_scalar_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.zero_scalar_accumulators();
                }
            }

            void zero_vector_accumulators() {
                domain.zero_vector_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.zero_vector_accumulators();
                }
            }

            void zero_matrix_accumulators() {
                domain.zero_matrix_accumulators();

                for (auto &p : boundary) {
                    auto &b = p.second;
                    b.zero_matrix_accumulators();
                }
            }

            template <class ForcingFunctionDescription>
            void add_forcing_function(const ForcingFunctionDescription &desc) {
                if (!domain.fe) {
                    assert(false);
                }

                auto assembler =
                    std::make_shared<utopia::intrepid2::Assemble<ForcingFunctionDescription>>(domain.fe, desc);

                domain.assemblers.push_back(assembler);
            }

            template <class ForcingFunctionDescription>
            void add_forcing_function_on_boundary(const std::string &boundary_name,
                                                  const ForcingFunctionDescription &desc) {
                auto &b = boundary[boundary_name];
                std::shared_ptr<FE> bfe;

                if (!b.fe) {
                    bfe = std::make_shared<FE>();
                    create_fe_on_boundary(*space, *bfe, boundary_name, 2);
                    b.fe = bfe;
                    b.name = boundary_name;
                } else {
                    bfe = b.fe;
                }

                auto assembler = std::make_shared<utopia::intrepid2::Assemble<ForcingFunctionDescription>>(bfe, desc);
                b.assemblers.push_back(assembler);
            }

            void compute_residuals_from_matrices() {
                // Compute linear residuals in one go, Only works if we exclusively have linear materials
                domain.ensure_vector_accumulator();

                if (domain.has_matrix()) {
                    residual(
                        domain.matrix_accumulator->data(), x_field->data(), domain.vector_accumulator->data(), mode);
                }
            }

            void update(const Vector &x) {
                ensure_field();

                assert(space);
                utopia::Field<FunctionSpace> in("x", space, make_ref(const_cast<Vector &>(x)));
                in.set_tensor_size(space->n_var());
                convert_field(in, *x_field);

                for (auto &a_ptr : domain.assemblers) {
                    a_ptr->update(x_field);
                }

                for (auto &p : boundary) {
                    auto &b = p.second;

                    for (auto &a_ptr : b.assemblers) {
                        a_ptr->update(x_field);
                    }
                }
            }

            bool assemble_matrix(Matrix &mat) {
                ensure_output(mat);
                ensure_matrix_accumulators();
                zero_matrix_accumulators();

                for (auto &a_ptr : domain.assemblers) {
                    if (a_ptr->is_matrix()) {
                        if (!a_ptr->assemble_matrix()) {
                            return false;
                        }
                    }
                }

                for (auto &p : boundary) {
                    auto &b = p.second;

                    for (auto &a_ptr : b.assemblers) {
                        if (a_ptr->is_matrix()) {
                            if (!a_ptr->assemble_matrix()) {
                                return false;
                            }
                        }
                    }
                }

                local_to_global_matrix_domain(mat);
                return true;
            }

            bool assemble_vector(Vector &vec) {
                ensure_output(vec);
                ensure_vector_accumulators();
                zero_vector_accumulators();

                for (auto &a_ptr : domain.assemblers) {
                    if (a_ptr->is_vector()) {
                        if (!a_ptr->assemble_vector()) {
                            return false;
                        }
                    }
                }

                for (auto &p : boundary) {
                    auto &b = p.second;

                    for (auto &a_ptr : b.assemblers) {
                        if (a_ptr->is_vector()) {
                            if (!a_ptr->assemble_vector()) {
                                return false;
                            }
                        }
                    }
                }

                local_to_global_vector_domain(vec);
                local_to_global_vector_boundary(vec);
                return true;
            }

            bool assemble_material(Matrix &mat, Vector &vec) {
                ensure_output(mat);
                ensure_output(vec);

                ensure_accumulators();
                zero_accumulators();

                for (auto &a_ptr : domain.assemblers) {
                    if (!a_ptr->assemble()) {
                        return false;
                    }
                }

                for (auto &p : boundary) {
                    auto &b = p.second;

                    for (auto &a_ptr : b.assemblers) {
                        if (!a_ptr->assemble()) {
                            return false;
                        }
                    }
                }

                local_to_global_matrix_domain(mat);
                local_to_global_vector_domain(vec);
                local_to_global_vector_boundary(vec);
                return true;
            }

            void local_to_global_matrix_domain(Matrix &mat) {
                if (domain.has_matrix()) {
                    local_to_global(*space, domain.matrix_accumulator->data(), mode, mat);
                }
            }

            void local_to_global_vector_domain(Vector &vec) {
                if (domain.has_vector()) {
                    local_to_global(*space, domain.vector_accumulator->data(), mode, vec);
                }
            }

            void local_to_global_vector_boundary(Vector &vec) {
                for (auto &p : boundary) {
                    auto &b = p.second;

                    if (b.has_vector()) {
                        side_local_to_global(*space, b.vector_accumulator->data(), mode, vec, b.name);
                    }

                    assert(!b.has_matrix() && "IMPLEMENT ME");
                }
            }

            void ensure_output(Matrix &mat) {
                if (!mat.empty() && mat.is_assembled()) {
                    mat *= 0.0;
                }
            }

            void ensure_output(Vector &vec) {
                if (!vec.empty()) {
                    vec *= 0.0;
                }
            }

            void ensure_field() {
                if (!x_field) {
                    assert(domain.fe);
                    x_field = std::make_shared<intrepid2::Field<Scalar>>(domain.fe);
                }
            }

            std::shared_ptr<FunctionSpace> space;

            // Volume (only one for now)
            DomainAssembler domain;

            // Boundary
            std::map<std::string, BoundaryAssembler> boundary;

            // Env and Utils
            std::shared_ptr<Environment<FunctionSpace>> env;
            AssemblerRegistry registry;

            std::shared_ptr<intrepid2::Field<Scalar>> x_field;

            AssemblyMode mode{ADD_MODE};
            bool is_linear_{true};
        };

        template <class FunctionSpace>
        AssemblyMode OmniAssembler<FunctionSpace>::mode() const {
            return impl_->mode;
        }
        template <class FunctionSpace>
        void OmniAssembler<FunctionSpace>::set_mode(AssemblyMode mode) {
            impl_->mode = mode;
        }

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
        bool OmniAssembler<FunctionSpace>::assemble(const Vector &x, Matrix &matrix, Vector &fun) {
            if (!impl_->domain.fe) {
                return false;
            }

            impl_->update(x);
            if (!impl_->assemble_material(matrix, fun)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::assemble(const Vector &x, Matrix &matrix) {
            if (!impl_->domain.fe) {
                return false;
            }

            impl_->update(x);
            if (!impl_->assemble_matrix(matrix)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::assemble(const Vector &x, Vector &vec) {
            if (!impl_->domain.fe) {
                return false;
            }

            impl_->update(x);
            if (!impl_->assemble_vector(vec)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::assemble(Matrix &matrix) {
            if (!impl_->domain.fe) {
                return false;
            }

            if (!impl_->assemble_matrix(matrix)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::assemble(Vector &vector) {
            if (!impl_->domain.fe) {
                return false;
            }

            if (!impl_->assemble_vector(vector)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace>
        void OmniAssembler<FunctionSpace>::read(Input &in) {
            // FIXME order must be guessed by discretization and material
            int quadrature_order = 2;
            in.get("quadrature_order", quadrature_order);
            impl_->domain.fe = std::make_shared<typename Impl::FE>();
            create_fe(*impl_->space, *impl_->domain.fe, quadrature_order);

            impl_->is_linear_ = true;

            in.get("material", [this](Input &node) {
                auto assembler = impl_->registry.make_assembler(impl_->domain.fe, node);

                if (assembler) {
                    if (!assembler->is_linear()) {
                        impl_->is_linear_ = false;
                    }

                    impl_->domain.assemblers.push_back(assembler);
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
                            impl_->add_forcing_function_on_boundary(name, ff);
                        }

                    } else {
                        if (forcing_function_type == "value") {
                            ForcingFunction<Scalar> ff;
                            ff.read(node);
                            ff.n_components = impl_->space->n_var();
                            impl_->add_forcing_function(ff);
                        }
                    }
                });
            });
        }

        template <class FunctionSpace>
        bool OmniAssembler<FunctionSpace>::is_linear() const {
            return impl_->is_linear_;
        }

    }  // namespace intrepid2
}  // namespace utopia
