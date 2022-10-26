#ifndef UTOPIA_MESH_TRANSFORM_HPP
#define UTOPIA_MESH_TRANSFORM_HPP

// basic
#include "utopia_Input.hpp"

// FE
#include "utopia_Field.hpp"

#include <algorithm>

namespace utopia {
    template <class FunctionSpace>
    class MeshTransform : public Configurable {
    public:
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;

        void read(Input &in) override {
            shift_.resize(3, 0);
            scale_factors_.resize(3, 1);

            in.get("xshift", shift_[0]);
            in.get("yshift", shift_[1]);
            in.get("zshift", shift_[2]);

            in.get("xscale", scale_factors_[0]);
            in.get("yscale", scale_factors_[1]);
            in.get("zscale", scale_factors_[2]);
        }

        void generate_displacement_field(FunctionSpace &space, Field<FunctionSpace> &displacement) {
            if (displacement.empty()) {
                space.create_field(displacement);
            }

            if (scale_factors_.empty()) {
                return;
            }

            int n_var = space.n_var();

            assert(n_var == space.mesh().spatial_dimension());

            auto d_view = view_device(displacement.data());

            Range r = range(displacement.data());
            SizeType r_begin = r.begin() / n_var;
            SizeType r_end = r.end() / n_var;

            const int n_factors = scale_factors_.size();

            int dim = std::min(n_factors, std::min(space.mesh().spatial_dimension(), n_var));

            auto fun = [=](const SizeType idx, const Scalar *point) {
                if (idx < r_begin || idx >= r_end) return;

                Scalar p3[3] = {0.0, 0.0, 0.0};
                Scalar transformed_p3[3] = {0.0, 0.0, 0.0};

                for (int d = 0; d < dim; ++d) {
                    p3[d] = point[d];
                }

                for (int d = 0; d < dim; ++d) {
                    transformed_p3[d] = (p3[d] + shift_[d]) * scale_factors_[d];
                }

                for (int d = 0; d < dim; ++d) {
                    d_view.set(idx * n_var + d, transformed_p3[d] - p3[d]);
                }
            };

            space.node_eval(fun);
        }

    private:
        std::vector<Scalar> shift_;
        std::vector<Scalar> scale_factors_;
    };
}  // namespace utopia

#endif  // UTOPIA_MESH_TRANSFORM_HPP
