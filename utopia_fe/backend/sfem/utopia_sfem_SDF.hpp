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
            using Vector = Traits<SDF>::Vector;

            SDF();
            ~SDF();
            void read(Input &in) override;
            void describe(std::ostream &os) const override;
            bool to_mesh(const Mesh &mesh, Vector &field);
            bool to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field);
            bool interpolate_to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field);

            bool project_to_mesh_with_sharp_features(const Mesh &mesh,
                                                     Vector &field,
                                                     Vector &grad_field,
                                                     Vector &weights);

            bool project_to_mesh(const Mesh &mesh, Vector &field, Vector &grad_field, Vector &weights);
            bool interpolate() const;
            bool has_weights() const;
            const Vector &weights() const;
            void clear();

            class Impl;
            std::unique_ptr<Impl> impl_;

        private:
            void read_from_file(const Path &path);
        };
    }  // namespace sfem
}  // namespace utopia

#endif  // UTOPIA_SFEM_SDF_HPP
