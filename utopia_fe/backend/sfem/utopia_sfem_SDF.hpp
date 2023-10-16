#ifndef UTOPIA_SFEM_SDF_HPP
#define UTOPIA_SFEM_SDF_HPP

#include "utopia_fe_base.hpp"

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_sfem_ForwardDeclarations.hpp"

namespace utopia {
    template <>
    class Traits<utopia::sfem::SDF> : public Traits<UVector> {
    public:
        using Super = Traits<UVector>;
    };

    namespace sfem {
        class SDF : public Configurable, public Describable {
        public:
            using Communicator = Traits<SDF>::Communicator;

            SDF();
            ~SDF();
            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace sfem
}  // namespace utopia

#endif  // UTOPIA_SFEM_SDF_HPP
