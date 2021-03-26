#include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_LinearElasticity.hpp"
#include "utopia_intrepid2_VectorLaplaceOperator.hpp"
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_intrepid2.hpp"

#include <functional>

namespace utopia {
    namespace stk {

        class OmniAssembler::Impl {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;

            template <class MaterialDescription>
            void init_assembler(const MaterialDescription &desc) {
                assemble = [this, desc](const Vector &x, Matrix &mat, Vector &rhs) -> bool {
                    utopia::intrepid2::Assemble<MaterialDescription> assembler(desc, fe);
                    assembler.init();
                    local_to_global(*space, assembler.element_matrices(), mat);

                    rhs = mat * x;
                    // rhs.set(0.0);
                    return true;
                };
            }

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

            return impl_->assemble(x, jacobian, fun);
        }

        void OmniAssembler::read(Input &in) {
            // FIXME order must be guessed by discretization and material
            int quadrature_order = 2;
            in.get("quadrature_order", quadrature_order);
            impl_->fe = std::make_shared<Impl::FE>();
            create_fe(*impl_->space, *impl_->fe, quadrature_order);

            std::string material_type = "";

            in.get("material", [&](Input &in) { in.get("type", material_type); });

            // FIXME create a registry and decentralize material registration
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

                return;
            }
        }

    }  // namespace stk
}  // namespace utopia
