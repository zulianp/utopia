#ifndef UTOPIA_MOONOLITH_CONTACT_HPP
#define UTOPIA_MOONOLITH_CONTACT_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_moonolith_ForwardDeclarations.hpp"
#include "utopia_moonolith_FunctionSpace.hpp"

#include <memory>
#include <unordered_set>

namespace utopia {
    namespace moonolith {

        class Contact final : public Configurable, public Describable {
        public:
            using Vector = Traits<FunctionSpace>::Vector;
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using IndexArray = Traits<FunctionSpace>::IndexArray;
            using Comm = Traits<FunctionSpace>::Communicator;

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            bool assemble(const FunctionSpace &space);
            void transform(const Matrix &in, Matrix &out);
            void transform(const Vector &in, Vector &out);
            void inverse_transform(const Vector &in, Vector &out);

            const Vector &gap() const;
            const Vector &is_contact() const;
            const Vector &normals() const;

            Vector &gap();
            Vector &is_contact();
            Vector &normals();

            std::shared_ptr<Matrix> mass_matrix();
            std::shared_ptr<Matrix> orthogonal_transformation();
            std::shared_ptr<Matrix> complete_transformation();

            Contact();
            ~Contact();

            class Params : public Configurable, public Describable {
            public:
                //     int variable_number{0};
                //     Scalar gap_negative_bound{-0.0001};
                //     Scalar gap_positive_bound{0.1};
                //     Scalar margin{0.};
                //     std::unordered_set<int> tags;
                //     bool invert_face_orientation{false};
                //     bool debug{false};
                //     bool snap_to_canonical_vectors{false};
                //     int skip_dir{-1};
                //     Scalar skip_dir_tol{0.2};
                //     bool verbose{false};

                void read(Input &in) override;
                void describe(std::ostream &os) const override;
            };

            void set_params(const Params &params);
            inline const Params &params() const { return *params_; }

            void set_banned_nodes(const std::shared_ptr<IndexArray> &banned_nodes);

            UTOPIA_NVCC_PRIVATE
            class Impl;
            class Output;

            template <int Dim>
            class ImplD;

            std::unique_ptr<Params> params_;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace moonolith
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_CONTACT_HPP
