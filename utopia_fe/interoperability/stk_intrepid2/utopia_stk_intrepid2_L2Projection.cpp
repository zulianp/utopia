
#include "utopia_stk_intrepid2_L2Projection.hpp"

#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_L2Projection.hpp"
#include "utopia_kokkos_Mass.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    namespace stk {

        class L2Projection::Impl {
        public:
            using Intrepi2FE_t = utopia::intrepid2::FE<Scalar>;
            using L2Projection2 = utopia::kokkos::L2Projection<FunctionSpace, Intrepi2FE_t, Intrepi2FE_t::DynRankView>;
            using Intrepid2Assembler = utopia::kokkos::FEAssembler<FunctionSpace, Intrepi2FE_t>;

            std::shared_ptr<Field> field;
            std::shared_ptr<Intrepid2Field> intrepid2_field;
            bool verbose{false};
            bool print_field{false};
        };

        void L2Projection::set_field(const std::shared_ptr<Field> &field) { impl_->field = field; }
        void L2Projection::set_field(const std::shared_ptr<Intrepid2Field> &field) { impl_->intrepid2_field = field; }

        L2Projection::L2Projection(const std::shared_ptr<FE> &fe) : Super(fe), impl_(utopia::make_unique<Impl>()) {}

        L2Projection::~L2Projection() = default;

        bool L2Projection::assemble_vector() {
            ensure_assembler();

            assert(this->assembler());
            return this->assembler()->assemble_vector();
        }

        void L2Projection::ensure_assembler() {
            if (!this->assembler()) {
                using Assembler = Impl::L2Projection2;
                auto assembler = std::make_shared<Assembler>(this->fe_ptr(), impl_->intrepid2_field->data());
                this->set_assembler(assembler);
            }
        }

        void L2Projection::read(Input &in) {
            Super::read(in);

            if (!impl_->intrepid2_field) {
                if (!impl_->field) {
                    auto env = this->environment();

                    if (!env) {
                        assert(false);
                        Utopia::Abort("Environment not defined in L2Projection");
                    }

                    std::string field;
                    in.require("field", field);

                    impl_->field = env->find_field(*this->space(), field);

                    if (!impl_->field) {
                        utopia::err() << "field not found in environment:";
                        env->describe(utopia::err().stream());
                    }
                }

                impl_->intrepid2_field = std::make_shared<Intrepid2Field>(this->fe_ptr());
                convert_field(*impl_->field, *impl_->intrepid2_field);
            }

            in.get("verbose", impl_->verbose);
            in.get("print_field", impl_->print_field);

            // const int spatial_dim = this->space()->mesh().spatial_dimension();

            if (impl_->print_field) {
                impl_->intrepid2_field->describe(utopia::out().stream());
            }

            ensure_assembler();
            this->assembler()->read(in);

            if (impl_->verbose) {
                utopia::out() << "-----------------------------\n";
                utopia::out() << "L2Projection\n";
                if (impl_->field) {
                    utopia::out() << "Field:\t" << impl_->field->name() << '\n';
                }

                utopia::out() << "print_field:\t" << impl_->print_field << '\n';
                utopia::out() << "-----------------------------\n";
            }
        }

    }  // namespace stk

}  // namespace utopia
