#ifndef UTOPIA_MOONOLITH_FUNCTION_SPACE_NEW_HPP
#define UTOPIA_MOONOLITH_FUNCTION_SPACE_NEW_HPP

#include "utopia_moonolith_Mesh.hpp"

namespace utopia {
    template <>
    class Traits<utopia::moonolith::FunctionSpace> : public Traits<utopia::moonolith::Mesh> {};

    namespace moonolith {

        class FunctionSpace : public Configurable, public Describable, public Traits<FunctionSpace> {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using Comm = Traits<FunctionSpace>::Communicator;

            FunctionSpace(const Comm &comm = Comm::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void init(const std::shared_ptr<Mesh> &mesh, const bool init_as_iso_paramatric = true);

            bool write(const Path &path, const Vector &x);
            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            std::shared_ptr<Mesh> mesh_ptr() const;
            const Mesh &mesh() const;
            Mesh &mesh();

            inline const Comm &comm() const { return mesh().comm(); }

            template <int Dim>
            std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<Scalar, Dim>>> raw_type() const;

            template <int Dim>
            void wrap(const std::shared_ptr<::moonolith::FunctionSpace<::moonolith::Mesh<Scalar, Dim>>> &space);

            SizeType n_dofs() const;
            SizeType n_local_dofs() const;

            // void create_matrix(Matrix &mat) const;
            // void create_vector(Vector &vec) const;

        private:
            class Impl;
            std::shared_ptr<Impl> impl_;
        };

    }  // namespace moonolith

}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_FUNCTION_SPACE_NEW_HPP
