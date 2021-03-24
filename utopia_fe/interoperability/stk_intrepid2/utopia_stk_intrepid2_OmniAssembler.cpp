#include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_LinearElasticity.hpp"
#include "utopia_intrepid2_VectorLaplaceOperator.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include <functional>

namespace utopia {
    namespace stk {

        class OmniAssembler::Impl {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;

            template <class MaterialDescription>
            void init_assembler(const MaterialDescription &desc) {}

            std::function<bool(const Vector &x, Matrix &, Vector &)> assemble;
            std::shared_ptr<stk::FunctionSpace> space;
            std::shared_ptr<FE> fe;
        };

        OmniAssembler::OmniAssembler(const std::shared_ptr<stk::FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {
            impl_->space = space;
        }

        OmniAssembler::~OmniAssembler() = default;

        bool OmniAssembler::assemble(const Vector &x, Matrix &jacobian, Vector &fun) {
            if (!impl_->fe) {
                return false;
            }

            // IMPLEMENT ME
            return true;
        }

        void OmniAssembler::read(Input &in) {
            std::string material_type = "";

            in.get("material", [&](Input &in) { in.get("type", material_type); });

            if (material_type == "LaplaceOperator") {
                LaplaceOperator<Scalar> material(1.0);
                in.get("material", material);
                impl_->init_assembler(material);
            } else if (material_type == "VectorLaplaceOperator1") {
                VectorLaplaceOperator<1, Scalar> material(1.0);
                in.get("material", material);
                impl_->init_assembler(material);
            } else if (material_type == "VectorLaplaceOperator2") {
                VectorLaplaceOperator<2, Scalar> material(1.0);
                in.get("material", material);
                impl_->init_assembler(material);
            } else if (material_type == "VectorLaplaceOperator3") {
                VectorLaplaceOperator<3, Scalar> material(1.0);
                in.get("material", material);
                impl_->init_assembler(material);
            } else if (material_type == "LinearElasticity1") {
                LinearElasticity<1, Scalar> material(1.0, 1.0);
                in.get("material", material);
                impl_->init_assembler(material);
            } else if (material_type == "LinearElasticity2") {
                LinearElasticity<2, Scalar> material(1.0, 1.0);
                in.get("material", material);
                impl_->init_assembler(material);
            } else if (material_type == "LinearElasticity3") {
                LinearElasticity<3, Scalar> material(1.0, 1.0);
                in.get("material", material);
                impl_->init_assembler(material);
            } else {
                if (material_type.empty()) {
                    utopia::err() << "[Error] Undefined material\n";
                } else {
                    utopia::err() << "[Error] Unsupported material " << material_type << '\n';
                }
            }
        }

    }  // namespace stk
}  // namespace utopia
