#ifndef UTOPIA_KOKKOS_MATERIAL_HPP
#define UTOPIA_KOKKOS_MATERIAL_HPP

#include "utopia_fe_Environment.hpp"
#include "utopia_kokkos_Discretization.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace utopia {
    namespace kokkos {
        template <class FunctionSpace, class FE>
        class FEKernelAssembler {
        public:
            void matrix_assembly_begin();
            void matrix_assembly_end();

            void vector_assembly_begin();
            void vector_assembly_end();

            void scalar_assembly_begin();
            void scalar_assembly_end();

            template <class Op>
            bool matrix_kernel();

            template <class Op>
            bool apply_kernel();

            template <class Op>
            bool vector_kernel();

            template <class Op>
            bool scalar_kernel();
        };

        template <class FunctionSpace, class FE, class Assembler>
        class Material : public Configurable {
        public:
            using Discretization = utopia::kokkos::Discretization<FunctionSpace, FE>;
            using Environment = utopia::Environment<FunctionSpace>;

            virtual ~Material() = default;

            virtual bool hessian(AssemblyMode mode) = 0;
            virtual bool gradient(AssemblyMode mode) = 0;
            virtual bool value(AssemblyMode mode) = 0;

            // Matrix free hessian application
            virtual bool apply(AssemblyMode mode) = 0;

            virtual void set_assembler(const std::shared_ptr<Assembler> &assembler) { assembler_ = assembler; }

            virtual void set_discretization(const std::shared_ptr<Assembler> &discretization) {
                discretization_ = discretization;
            }

            virtual void set_environment(const std::shared_ptr<Assembler> &environment) { environment_ = environment; }

            std::vector<std::shared_ptr<Field>> field(const std::string &name,
                                                      const typename Discretization::Part part = Discretization::all) {
                auto it = fields_.find(name);
                if (it != fields_.end()) return it->second;

                if (!discretization_) {
                    assert(false);
                    Utopia::Abort();
                }

                auto space = discretization_->space();

                if (!space_) {
                    assert(false);
                    Utopia::Abort();
                }

                auto f = environment_->find_field(name, space);

                auto &fef = fields_[name];
                discretization_->convert_field(f, fef, part);
                return fef;
            }

            void clear_cache() { fields_.clear(); }

            void read(Input &) override {
                // TODO
            }

        private:
            std::shared_ptr<Assembler> assembler_;
            std::shared_ptr<Discretization> discretization_;
            std::shared_ptr<Environment> environment_;
            std::map<std::string, std::vector<std::shared_ptr<Field>>> fields_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MATERIAL_HPP
