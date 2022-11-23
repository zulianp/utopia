#ifndef UTOPIA_KOKKOS_MATERIAL_HPP
#define UTOPIA_KOKKOS_MATERIAL_HPP

#include "utopia_Function.hpp"

#include "utopia_fe_Environment.hpp"
#include "utopia_kokkos_Discretization.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE, class Assembler>
        class Material
            : public utopia::Function<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector>,
              public Operator<typename Traits<FunctionSpace>::Vector> {
        public:
            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using Communicator = typename Traits<FunctionSpace>::Communicator;

            using Super = utopia::Function<Matrix, Vector>;

            using Discretization = utopia::kokkos::Discretization<FunctionSpace, FE>;
            using Environment = utopia::Environment<FunctionSpace>;
            using Field = utopia::kokkos::Field<FE>;

            virtual ~Material() = default;

            virtual int n_vars() const = 0;
            virtual std::string name() const = 0;

            virtual bool has_hessian() const = 0;
            virtual bool has_gradient() const = 0;
            virtual bool has_value() const = 0;
            virtual bool is_operator() const = 0;
            virtual bool is_linear() const { return false; }

            bool is_hessian_constant() const override { return is_linear(); }

            Size size() const override {
                auto n = assembler()->discretization()->space()->n_dofs();
                return {n, n};
            }

            Size local_size() const override {
                auto nl = assembler()->discretization()->space()->n_local_dofs();
                return {nl, nl};
            }

            const Communicator &comm() const override { return assembler()->discretization()->space()->comm(); }

            virtual bool apply(const Vector &x, Vector &y) const {
                if (!is_operator()) {
                    assert(false);
                    return false;
                }

                assembler()->vector_assembly_begin(y, mode_);

                assembler()->update_input(x);

                auto sol = assembler()->current_solution();

                bool ok = const_cast<Material *>(this)->apply_assemble(*sol, mode_);
                if (!ok) {
                    return false;
                }

                assembler()->vector_assembly_end(y, mode_);
                return false;
            }

            virtual bool value(const Vector &x, Scalar &value) const {
                if (!has_value()) {
                    assert(false);
                    return false;
                }

                assembler()->scalar_assembly_begin(value, mode_);

                assembler()->update_input(x);

                bool ok = const_cast<Material *>(this)->value_assemble(mode_);
                if (!ok) {
                    return false;
                }

                assembler()->scalar_assembly_end(value, mode_);
                return true;
            }

            virtual bool gradient(const Vector &x, Vector &g) const {
                if (!has_gradient() && is_linear() && is_operator()) {
                } else if (!has_gradient()) {
                    assert(false);
                    return false;
                }

                assembler()->vector_assembly_begin(g, mode_);

                assembler()->update_input(x);

                if (has_gradient()) {
                    bool ok = const_cast<Material *>(this)->gradient_assemble(mode_);
                    if (!ok) {
                        return false;
                    }
                }

                assembler()->vector_assembly_end(g, mode_);
                return true;
            }

            virtual bool hessian(const Vector &x, Matrix &H) const {
                if (!has_hessian()) {
                    assert(false);
                    return false;
                }

                assembler()->matrix_assembly_begin(H, mode_);

                assembler()->update_input(x);

                bool ok = const_cast<Material *>(this)->hessian_assemble(mode_);
                if (!ok) {
                    return false;
                }

                assembler()->matrix_assembly_end(H, mode_);
                return true;
            }

            virtual bool hessian_assemble(AssemblyMode mode) = 0;
            virtual bool gradient_assemble(AssemblyMode) { return false; }
            virtual bool value_assemble(AssemblyMode mode) = 0;

            // Matrix free hessian application
            virtual bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) = 0;

            virtual void set_assembler(const std::shared_ptr<Assembler> &assembler) { assembler_ = assembler; }

            inline std::shared_ptr<Assembler> assembler() const {
                assert(assembler_);
                return assembler_;
            }

            virtual void set_environment(const std::shared_ptr<Environment> &environment) {
                environment_ = environment;
            }

            std::vector<std::shared_ptr<Field>> field(const std::string &name,
                                                      const typename Discretization::Part part = Discretization::all) {
                auto it = fields_.find(name);
                if (it != fields_.end()) return it->second;

                if (!assembler_) {
                    assert(false);
                    Utopia::Abort();
                }

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

                auto f = environment_->find_field(name, space);

                auto &fef = fields_[name];
                discretization->convert_field(f, fef, part);
                return fef;
            }

            void clear_cache() { fields_.clear(); }

            void read(Input &) override {
                // TODO
            }

        private:
            AssemblyMode mode_{OVERWRITE_MODE};
            std::shared_ptr<Assembler> assembler_;
            std::shared_ptr<Environment> environment_;
            std::map<std::string, std::vector<std::shared_ptr<Field>>> fields_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MATERIAL_HPP
