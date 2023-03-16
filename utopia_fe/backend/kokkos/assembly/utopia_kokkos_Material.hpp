#ifndef UTOPIA_KOKKOS_MATERIAL_HPP
#define UTOPIA_KOKKOS_MATERIAL_HPP

#include "utopia_Material.hpp"

#include "utopia_fe_Environment.hpp"
#include "utopia_kokkos_Discretization.hpp"

#include "utopia_kokkos_Gradient.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE, class Assembler>
        class Material : public utopia::AbstractMaterial<FunctionSpace> {
        public:
            using Super = utopia::AbstractMaterial<FunctionSpace>;

            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using Communicator = typename Traits<FunctionSpace>::Communicator;

            using Discretization = utopia::kokkos::Discretization<FunctionSpace, FE>;
            using Environment = utopia::Environment<FunctionSpace>;
            using Field = utopia::kokkos::Field<FE>;
            using FieldPtr = std::shared_ptr<utopia::kokkos::Field<FE>>;

            virtual ~Material() = default;

            virtual int n_vars() const = 0;

            virtual int order() const { return 2; }

            Size size() const override {
                auto n = assembler()->discretization()->space()->n_dofs();
                return {n, n};
            }

            Size local_size() const override {
                auto nl = assembler()->discretization()->space()->n_local_dofs();
                return {nl, nl};
            }

            const Communicator &comm() const override { return assembler()->discretization()->space()->comm(); }

            bool apply(const Vector &x, Vector &y) const override {
                if (!this->is_operator()) {
                    assert(false);
                    return false;
                }

                UTOPIA_TRACE_REGION_BEGIN("utopia::kokkos::Material::apply");

                assembler()->vector_assembly_begin(y, mode_);

                assembler()->update_input(x);
                solution_updated(assembler()->current_solution());
                auto sol = assembler()->current_solution();

                bool ok = const_cast<Material *>(this)->apply_assemble(*sol, mode_);
                assembler()->vector_assembly_end(y, mode_);

                if (must_apply_constraints_) {
                    assembler()->discretization()->space()->copy_at_constrained_nodes(x, y);
                }

                UTOPIA_TRACE_REGION_END("utopia::kokkos::Material::apply");
                return ok;
            }

            bool value(const Vector &x, Scalar &value) const override {
                if (!this->has_value()) {
                    assert(false);
                    return false;
                }

                UTOPIA_TRACE_REGION_BEGIN("utopia::kokkos::Material::value");

                // if (OVERWRITE_MODE == mode_)
                // value = 0;

                assembler()->scalar_assembly_begin(value, mode_);

                assembler()->update_input(x);
                solution_updated(assembler()->current_solution());

                bool ok = const_cast<Material *>(this)->value_assemble(mode_);
                assembler()->scalar_assembly_end(value, mode_);

                UTOPIA_TRACE_REGION_END("utopia::kokkos::Material::value");
                return ok;
            }

            bool gradient(const Vector &x, Vector &g) const override {
                if (!this->has_gradient() && this->is_linear() && this->is_operator()) {
                } else if (!this->has_gradient()) {
                    assert(false);
                    return false;
                }

                UTOPIA_TRACE_REGION_BEGIN("utopia::kokkos::Material::gradient");

                assembler()->vector_assembly_begin(g, mode_);

                assembler()->update_input(x);
                solution_updated(assembler()->current_solution());

                bool ok = false;
                if (this->has_gradient()) {
                    ok = const_cast<Material *>(this)->gradient_assemble(mode_);
                } else {
                    ok = const_cast<Material *>(this)->apply_assemble(*assembler()->current_solution(), mode_);
                }

                assembler()->vector_assembly_end(g, mode_);

                if (must_apply_constraints_) {
                    assembler()->discretization()->space()->apply_zero_constraints(g);
                }

                UTOPIA_TRACE_REGION_END("utopia::kokkos::Material::gradient");
                return ok;
            }

            virtual bool hessian(const Vector &x, Matrix &H) const {
                if (!this->has_hessian()) {
                    assert(false);
                    return false;
                }

                UTOPIA_TRACE_REGION_BEGIN("utopia::kokkos::Material::hessian");

                assembler()->matrix_assembly_begin(H, mode_);

                assembler()->update_input(x);
                solution_updated(assembler()->current_solution());

                bool ok = const_cast<Material *>(this)->hessian_assemble(mode_);

                assembler()->matrix_assembly_end(H, mode_);

                if (must_apply_constraints_) {
                    assembler()->discretization()->space()->apply_constraints(H);
                }

                UTOPIA_TRACE_REGION_END("utopia::kokkos::Material::hessian");
                return ok;
            }

            virtual bool hessian_assemble(AssemblyMode mode) = 0;
            virtual bool gradient_assemble(AssemblyMode) { return false; }
            virtual bool value_assemble(AssemblyMode mode) = 0;

            virtual bool solution_updated(const FieldPtr &) const { return true; }

            // Matrix free hessian application
            virtual bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) = 0;

            virtual void set_assembler(const std::shared_ptr<Assembler> &assembler) {
                assembler_ = assembler;
                assembler_->set_n_vars(this->n_vars());
            }

            inline std::shared_ptr<Assembler> assembler() const {
                assert(assembler_);
                return assembler_;
            }

            virtual void set_environment(const std::shared_ptr<Environment> &environment) {
                environment_ = environment;
            }

            std::vector<std::shared_ptr<Field>> field(
                const std::string &name,
                const typename Discretization::Part part = Discretization::all()) {
                auto it = fields_.find(name);
                if (it != fields_.end()) return it->second;

                if (!assembler_) {
                    assert(false);
                    Utopia::Abort();
                }

                UTOPIA_TRACE_REGION_BEGIN("utopia::kokkos::Material::field");

                auto discretization = assembler()->discretization();

                if (!discretization) {
                    assert(false);
                    Utopia::Abort();
                }

                auto space = discretization->space();

                if (!space) {
                    assert(false);
                    Utopia::Abort();
                }

                auto f = environment_->find_field(*space, name);

                if (!f) {
                    return {};
                }

                auto &fef = fields_[name];
                discretization->convert_field(*f, fef, part);

                UTOPIA_TRACE_REGION_END("utopia::kokkos::Material::field");
                return fef;
            }

            void clear_cache() { fields_.clear(); }

            void read(Input &) override {
                // TODO
            }

            // virtual void initialize(const std::shared_ptr<FunctionSpace> &space) = 0;

            // virtual void initialize_boundary(const std::shared_ptr<FunctionSpace> &space,
            //                                  const std::string &boundary_name) {
            //     Utopia::Abort();
            // }

            // virtual void initialize_block(const std::shared_ptr<FunctionSpace> &space, const std::string &block_name)
            // {
            //     Utopia::Abort();
            // }

        private:
            AssemblyMode mode_{OVERWRITE_MODE};
            std::shared_ptr<Assembler> assembler_;
            std::shared_ptr<Environment> environment_;
            std::map<std::string, std::vector<std::shared_ptr<Field>>> fields_;
            bool must_apply_constraints_{true};
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MATERIAL_HPP
