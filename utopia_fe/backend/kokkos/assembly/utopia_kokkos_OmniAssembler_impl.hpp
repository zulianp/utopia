#include "utopia_kokkos_OmniAssembler.hpp"

// utopia includes
#include "utopia_CreateFE.hpp"
#include "utopia_LocalToGlobal.hpp"

// utopia/kokkos
#include "utopia_kokkos_AutoHyperElasticity.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_ForcingFunction.hpp"
#include "utopia_kokkos_IncrementalForcingFunction.hpp"
#include "utopia_kokkos_IsotropicPhaseFieldForBrittleFractures.hpp"
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_LinearElasticity.hpp"
#include "utopia_kokkos_Mass.hpp"
#include "utopia_kokkos_NeoHookean.hpp"
#include "utopia_kokkos_Residual.hpp"
#include "utopia_kokkos_Transport.hpp"
#include "utopia_kokkos_VectorLaplaceOperator.hpp"
#include "utopia_kokkos_WeakLinearThermoElasticity.hpp"

// utopia/kokkos includes
#include "utopia_kokkos_FE.hpp"

#include <functional>

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE>
        class AssemblerRegistry {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using FE_t = FE;
            using FEAssembler_t = utopia::kokkos::FEAssembler<FE_t>;
            using FEAssemblerPtr_t = std::shared_ptr<FEAssembler_t>;

            FEAssemblerPtr_t make_assembler(const std::shared_ptr<FE_t> &fe, Input &in) const {
                std::string type;
                in.get("type", type);

                if (has_variants(type)) {
                    type = create_variant_name(type, fe->spatial_dimension());
                }

                auto it = assemblers_.find(type);
                if (it == assemblers_.end()) {
                    return nullptr;
                } else {
                    return it->second(fe, in);
                }
            }

            template <typename MaterialDescription>
            void register_assembler(const std::string &type) {
                assemblers_[type] = [](const std::shared_ptr<FE_t> &fe, Input &in) -> FEAssemblerPtr_t {
                    typename MaterialDescription::Params mat_desc;
                    mat_desc.read(in);
                    return std::make_shared<MaterialDescription>(fe, mat_desc);
                };
            }

            template <typename MaterialDescription>
            void register_assembler_variant(const std::string &type, const int spatial_dimension) {
                assemblers_[create_variant_name(type, spatial_dimension)] = [](const std::shared_ptr<FE_t> &fe,
                                                                               Input &in) -> FEAssemblerPtr_t {
                    typename MaterialDescription::Params mat_desc;
                    mat_desc.read(in);
                    return std::make_shared<MaterialDescription>(fe, mat_desc);
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
                register_assembler<utopia::kokkos::Mass<FE_t>>("Mass");
                register_assembler<utopia::kokkos::LaplaceOperator<FE_t>>("LaplaceOperator");
                register_assembler<utopia::kokkos::ForcingFunction<FE_t>>("ForcingFunction");
                register_assembler<utopia::kokkos::IncrementalForcingFunction<FE_t>>("IncrementalForcingFunction");
                register_assembler<utopia::kokkos::NeoHookean<FE_t>>("NeoHookean");
                register_assembler<utopia::kokkos::AutoHyperElasticity<FE_t>>("AutoHyperElasticity");

                register_assembler_variant<utopia::kokkos::VectorLaplaceOperator<FE_t, 1, Scalar_t>>(
                    "VectorLaplaceOperator", 1);
                register_assembler_variant<utopia::kokkos::VectorLaplaceOperator<FE_t, 2, Scalar_t>>(
                    "VectorLaplaceOperator", 2);
                register_assembler_variant<utopia::kokkos::VectorLaplaceOperator<FE_t, 3, Scalar_t>>(
                    "VectorLaplaceOperator", 3);

                register_assembler_variant<utopia::kokkos::LinearElasticity<FE_t, 1, Scalar_t>>("LinearElasticity", 1);
                register_assembler_variant<utopia::kokkos::LinearElasticity<FE_t, 2, Scalar_t>>("LinearElasticity", 2);
                register_assembler_variant<utopia::kokkos::LinearElasticity<FE_t, 3, Scalar_t>>("LinearElasticity", 3);

                register_assembler_variant<utopia::kokkos::WeakLinearThermoElasticity<FE_t, 1, Scalar_t>>(
                    "WeakLinearThermoElasticity", 1);
                register_assembler_variant<utopia::kokkos::WeakLinearThermoElasticity<FE_t, 2, Scalar_t>>(
                    "WeakLinearThermoElasticity", 2);
                register_assembler_variant<utopia::kokkos::WeakLinearThermoElasticity<FE_t, 3, Scalar_t>>(
                    "WeakLinearThermoElasticity", 3);

                // register_assembler_variant<utopia::kokkos::IsotropicPhaseFieldForBrittleFractures<FE_t, 1,
                // Scalar_t>>(
                //     "IsotropicPhaseFieldForBrittleFractures", 1);
                register_assembler_variant<utopia::kokkos::IsotropicPhaseFieldForBrittleFractures<FE_t, 2>>(
                    "IsotropicPhaseFieldForBrittleFractures", 2);
                register_assembler_variant<utopia::kokkos::IsotropicPhaseFieldForBrittleFractures<FE_t, 3>>(
                    "IsotropicPhaseFieldForBrittleFractures", 3);
            }
        };

        template <class FunctionSpace, class FE_>
        class OmniAssembler<FunctionSpace, FE_>::Impl {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using FE = FE_;
            using Field = utopia::kokkos::Field<FE>;
            using AssemblerRegistry = utopia::kokkos::AssemblerRegistry<FunctionSpace, FE>;
            using FEAssembler_t = utopia::kokkos::FEAssembler<FE>;
            using FEAssemblerPtr_t = std::shared_ptr<FEAssembler_t>;
            using MatrixAccumulator_t = typename FEAssembler_t::MatrixAccumulator;
            using VectorAccumulator_t = typename FEAssembler_t::VectorAccumulator;
            using ScalarAccumulator_t = typename FEAssembler_t::ScalarAccumulator;

            class PartAssembler {
            public:
                virtual ~PartAssembler() {}

                std::string name;
                std::shared_ptr<FE> fe;
                std::vector<FEAssemblerPtr_t> assemblers;

                std::shared_ptr<MatrixAccumulator_t> matrix_accumulator;
                std::shared_ptr<VectorAccumulator_t> vector_accumulator;
                std::shared_ptr<ScalarAccumulator_t> scalar_accumulator;

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
            void add_forcing_function(const typename ForcingFunctionDescription::Params &desc) {
                if (!domain.fe) {
                    assert(false);
                }

                auto assembler = std::make_shared<ForcingFunctionDescription>(domain.fe, desc);

                domain.assemblers.push_back(assembler);
            }

            template <class ForcingFunctionDescription>
            void add_forcing_function_on_boundary(const std::string &boundary_name,
                                                  const int quadrature_order,
                                                  const typename ForcingFunctionDescription::Params &desc) {
                auto &b = boundary[boundary_name];
                std::shared_ptr<FE> bfe;

                if (!b.fe) {
                    bfe = std::make_shared<FE>();
                    create_fe_on_boundary(*space, *bfe, boundary_name, quadrature_order);
                    b.fe = bfe;
                    b.name = boundary_name;
                } else {
                    bfe = b.fe;
                }

                auto assembler = std::make_shared<ForcingFunctionDescription>(bfe, desc);
                b.assemblers.push_back(assembler);
            }

            void compute_residuals_from_matrices() {
                // Compute linear residuals in one go, Only works if we exclusively have linear materials
                domain.ensure_vector_accumulator();

                if (domain.has_matrix()) {
                    utopia::kokkos::residual(
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
                            assert(false);
                            return false;
                        }
                    }
                }

                for (auto &p : boundary) {
                    auto &b = p.second;

                    for (auto &a_ptr : b.assemblers) {
                        if (a_ptr->is_matrix()) {
                            if (!a_ptr->assemble_matrix()) {
                                assert(false);
                                return false;
                            }
                        }
                    }
                }

                local_to_global_matrix_domain(mat);
                return true;
            }

            bool apply(Vector &vec) {
                ensure_output(vec);
                ensure_vector_accumulators();
                zero_vector_accumulators();

                auto data = domain.vector_accumulator->data();

                for (auto &a_ptr : domain.assemblers) {
                    if (a_ptr->is_operator()) {
                        if (!a_ptr->apply(x_field->data(), data)) {
                            assert(false);
                            return false;
                        }
                    }
                }

                for (auto &p : boundary) {
                    auto &b = p.second;

                    for (auto &a_ptr : b.assemblers) {
                        if (a_ptr->is_operator()) {
                            assert(false);  // Not supported yet!
                            // if (!a_ptr->apply(*x_field, data)) {
                            //     assert(false);
                            //     return false;
                            // }
                        }
                    }
                }

                local_to_global_vector_domain(vec);
                // local_to_global_vector_boundary(vec);
                return true;
            }

            bool assemble_vector(Vector &vec) {
                ensure_output(vec);
                ensure_vector_accumulators();
                zero_vector_accumulators();

                for (auto &a_ptr : domain.assemblers) {
                    if (a_ptr->is_vector()) {
                        if (!a_ptr->assemble_vector()) {
                            assert(false);
                            return false;
                        }
                    }
                }

                for (auto &p : boundary) {
                    auto &b = p.second;

                    for (auto &a_ptr : b.assemblers) {
                        if (a_ptr->is_vector()) {
                            if (!a_ptr->assemble_vector()) {
                                assert(false);
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
                            assert(false);
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

                    if (ensure_scalar_matrix) {
                        if (mat.is_block()) {
                            mat.convert_to_scalar_matrix();
                        }
                    }
                }
            }

            void local_to_global_vector_domain(Vector &vec) {
                if (domain.has_vector()) {
                    local_to_global(*space, domain.vector_accumulator->data(), mode, vec);

                    // static bool first = true;

                    // if (first) {
                    //     space->write("boundary.e", vec);
                    //     first = false;
                    // }
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

                // Scalar sum_vec = sum(vec);
                // utopia::out() << "sum_vec: " << sum_vec << '\n';
                // disp(vec);

                // static bool first = true;

                // if (first) {
                //     space->write("boundary.e", vec);
                //     first = false;
                // }
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
                    x_field = std::make_shared<Field>(domain.fe);
                }
            }

            std::shared_ptr<FunctionSpace> space;

            // Volume (only one for now)
            DomainAssembler domain;

            // Boundary
            std::map<std::string, BoundaryAssembler> boundary;

            // Env and Utils
            std::shared_ptr<Environment> env;
            AssemblerRegistry registry;

            std::shared_ptr<Field> x_field;

            AssemblyMode mode{ADD_MODE};
            bool is_linear_{true};
            bool ensure_scalar_matrix{true};
            bool fail_if_unregistered{true};
            bool debug{false};
        };

        template <class FunctionSpace, class FE>
        AssemblyMode OmniAssembler<FunctionSpace, FE>::mode() const {
            return impl_->mode;
        }
        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::set_mode(AssemblyMode mode) {
            impl_->mode = mode;
        }

        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::set_environment(const std::shared_ptr<Environment> &env) {
            impl_->env = env;
        }

        template <class FunctionSpace, class FE>
        OmniAssembler<FunctionSpace, FE>::OmniAssembler(const std::shared_ptr<FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {
            impl_->space = space;
        }

        template <class FunctionSpace, class FE>
        OmniAssembler<FunctionSpace, FE>::~OmniAssembler() = default;

        template <class FunctionSpace, class FE>
        bool OmniAssembler<FunctionSpace, FE>::assemble(const Vector &x, Matrix &matrix, Vector &fun) {
            if (!impl_->domain.fe) {
                return false;
            }

            impl_->update(x);
            if (!impl_->assemble_material(matrix, fun)) {
                return false;
            }

            if (impl_->debug && fun.has_nan_or_inf()) {
                impl_->space->write("NaN.e", x);
                this->~OmniAssembler();
                assert(false);
                utopia::Utopia::Abort("Detected NaN in material gradient!");
            }

            return true;
        }

        template <class FunctionSpace, class FE>
        bool OmniAssembler<FunctionSpace, FE>::assemble(const Vector &x, Matrix &matrix) {
            if (!impl_->domain.fe) {
                return false;
            }

            impl_->update(x);
            if (!impl_->assemble_matrix(matrix)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace, class FE>
        bool OmniAssembler<FunctionSpace, FE>::assemble(const Vector &x, Vector &vec) {
            if (!impl_->domain.fe) {
                return false;
            }

            impl_->update(x);
            if (!impl_->assemble_vector(vec)) {
                return false;
            }

            if (impl_->debug && vec.has_nan_or_inf()) {
                impl_->space->write("NaN.e", x);
                this->~OmniAssembler();
                assert(false);
                utopia::Utopia::Abort("Detected NaN in material gradient!");
            }

            return true;
        }

        template <class FunctionSpace, class FE>
        bool OmniAssembler<FunctionSpace, FE>::apply(const Vector &x, Vector &hessian_times_x) {
            if (!impl_->domain.fe) {
                return false;
            }

            impl_->update(x);
            if (!impl_->apply(hessian_times_x)) {
                return false;
            }

            if (impl_->debug && hessian_times_x.has_nan_or_inf()) {
                impl_->space->write("NaN.e", x);
                this->~OmniAssembler();
                utopia::Utopia::Abort("Detected NaN in application of operator!");
            }

            return true;
        }

        template <class FunctionSpace, class FE>
        bool OmniAssembler<FunctionSpace, FE>::assemble(Matrix &matrix) {
            if (!impl_->domain.fe) {
                return false;
            }

            if (!impl_->assemble_matrix(matrix)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace, class FE>
        bool OmniAssembler<FunctionSpace, FE>::assemble(Vector &vector) {
            if (!impl_->domain.fe) {
                return false;
            }

            if (!impl_->assemble_vector(vector)) {
                return false;
            }

            return true;
        }

        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::add_domain_assembler(const Intrepid2FEAssemblerPtr &assembler) {
            impl_->domain.assemblers.push_back(assembler);
            if (!assembler->is_linear()) {
                impl_->is_linear_ = false;
            }
        }

        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::fail_if_unregistered(const bool val) {
            impl_->fail_if_unregistered = val;
        }

        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::set_domain_fe(const std::shared_ptr<FE> &fe) {
            impl_->domain.fe = fe;
        }

        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::read(Input &in) {
            using ForcingFunction_t = utopia::kokkos::ForcingFunction<typename Impl::FE>;
            using IncrementalForcingFunction_t = utopia::kokkos::IncrementalForcingFunction<typename Impl::FE>;

            if (!impl_->domain.fe) {
                // FIXME order must be guessed by discretization and material
                int quadrature_order = 2;
                in.get("quadrature_order", quadrature_order);
                impl_->domain.fe = std::make_shared<typename Impl::FE>();
                create_fe(*impl_->space, *impl_->domain.fe, quadrature_order);
            }

            impl_->is_linear_ = true;

            in.get("ensure_scalar_matrix", impl_->ensure_scalar_matrix);
            in.get("debug", impl_->debug);

            in.get("material", [this](Input &node) {
                auto assembler = impl_->registry.make_assembler(impl_->domain.fe, node);
                if (assembler) {
                    add_domain_assembler(assembler);
                } else if (impl_->fail_if_unregistered) {
                    assert(false && "Should not come here");
                    Utopia::Abort("Could not find material!");
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

                        int quadrature_order = 2;
                        node.get("quadrature_order", quadrature_order);

                        if (forcing_function_type == "value") {
                            typename ForcingFunction_t::Params ff;
                            ff.read(node);
                            ff.n_components = impl_->space->n_var();
                            impl_->template add_forcing_function_on_boundary<ForcingFunction_t>(
                                name, quadrature_order, ff);
                        } else if (forcing_function_type == "IncrementalForcingFunction") {
                            typename IncrementalForcingFunction_t::Params ff;
                            ff.read(node);
                            ff.n_components = impl_->space->n_var();
                            impl_->template add_forcing_function_on_boundary<IncrementalForcingFunction_t>(
                                name, quadrature_order, ff);
                        }

                    } else {
                        if (forcing_function_type == "value") {
                            typename ForcingFunction_t::Params ff;
                            ff.read(node);
                            ff.n_components = impl_->space->n_var();
                            impl_->template add_forcing_function<ForcingFunction_t>(ff);
                        } else if (forcing_function_type == "IncrementalForcingFunction") {
                            typename IncrementalForcingFunction_t::Params ff;
                            ff.read(node);
                            ff.n_components = impl_->space->n_var();
                            impl_->template add_forcing_function<IncrementalForcingFunction_t>(ff);
                        }
                    }
                });
            });
        }

        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::set_time(const std::shared_ptr<SimulationTime> &time) {
            for (auto &a : impl_->domain.assemblers) {
                a->set_time(time);
            }

            for (auto &p : impl_->boundary) {
                auto &b = p.second;

                for (auto &a : b.assemblers) {
                    a->set_time(time);
                }
            }
        }

        template <class FunctionSpace, class FE>
        std::string OmniAssembler<FunctionSpace, FE>::name() const {
            return "utopia::kokkos::OmniAssembler";
        }

        template <class FunctionSpace, class FE>
        bool OmniAssembler<FunctionSpace, FE>::is_linear() const {
            return impl_->is_linear_;
        }

        template <class FunctionSpace, class FE>
        std::shared_ptr<Environment<FunctionSpace>> OmniAssembler<FunctionSpace, FE>::environment() const {
            return impl_->env;
        }

        template <class FunctionSpace, class FE>
        void OmniAssembler<FunctionSpace, FE>::set_space(const std::shared_ptr<FunctionSpace> &space) {
            impl_->space = space;
        }

        template <class FunctionSpace, class FE>
        std::shared_ptr<FunctionSpace> OmniAssembler<FunctionSpace, FE>::space() const {
            return impl_->space;
        }

    }  // namespace kokkos
}  // namespace utopia
