
#include "utopia_stk_intrepid2_L2Projection.hpp"

#include "utopia_intrepid2_Gradient.hpp"
#include "utopia_intrepid2_L2Projection.hpp"
#include "utopia_intrepid2_Mass.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    namespace stk {

        class L2Projection::Impl {
        public:
            using L2Projection2 = utopia::L2Projection<Intrepid2FE::DynRankView>;
            using Intrepid2Assembler = utopia::intrepid2::FEAssembler<Scalar>;

            std::shared_ptr<Field> field;
            bool verbose{false};
            bool print_field{false};
        };

        void L2Projection::set_field(const std::shared_ptr<Field> &field) { impl_->field = field; }

        L2Projection::L2Projection(const std::shared_ptr<FE> &fe) : Super(fe), impl_(utopia::make_unique<Impl>()) {}

        L2Projection::~L2Projection() = default;

        void L2Projection::read(Input &in) {
            Super::read(in);

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

            in.get("verbose", impl_->verbose);
            in.get("print_field", impl_->print_field);

            const int spatial_dim = this->space()->mesh().spatial_dimension();

            Intrepid2Field intrepid2_field(this->fe_ptr());
            convert_field(*impl_->field, intrepid2_field);

            if (impl_->print_field) {
                intrepid2_field.describe(utopia::out().stream());
            }

            using Assembler = utopia::intrepid2::Assemble<Impl::L2Projection2>;
            auto assembler = std::make_shared<Assembler>(this->fe_ptr(), intrepid2_field.data());
            assembler->read(in);
            this->set_assembler(assembler);

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
