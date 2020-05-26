#ifndef UTOPIA_LIBMESH_INTERPOLATE_HPP
#define UTOPIA_LIBMESH_INTERPOLATE_HPP

#include "utopia_Interpolate.hpp"
#include "utopia_ProductFunctionSpace.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FEBackend.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include <vector>

namespace utopia {

    template <class FunSpace>
    using ProductTrialFunction = utopia::TrialFunction<ProductFunctionSpace<FunSpace>>;

    template <>
    class Interpolate<UVector, ProductTrialFunction<LibMeshFunctionSpace>>
        : public Expression<Interpolate<UVector, ProductTrialFunction<LibMeshFunctionSpace>>> {
    public:
        using Fun = utopia::ProductTrialFunction<LibMeshFunctionSpace>;
        using FunSpaceT = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;
        using FETraitsT = utopia::Traits<LibMeshFunctionSpace>;
        using Scalar = Traits<UVector>::Scalar;
        using ElementVector = FETraitsT::Vector;
        using ElementMatrix = FETraitsT::Matrix;
        using JacobianType = FETraitsT::JacobianType;

        using DenseVectorT = libMesh::DenseVector<libMesh::Real>;
        using FEBackend = utopia::FEBackend<LIBMESH_TAG>;

        class Data {
        public:
            Data(const UVector &coeff, Fun &fun) : coeff_(coeff), fun_(fun) {}

            const UVector &coeff_;
            UTOPIA_STORE(Fun) fun_;
            std::vector<std::unique_ptr<libMesh::FEBase>> fe_;
            std::vector<double> element_coeffs_;

            std::vector<DenseVectorT> fun_values_;
            std::vector<ElementMatrix> grad_values_;

            void fun_aux(std::vector<std::vector<DenseVectorT>> &ret) const {
                assert(!fe_.empty());

                auto &space = *fun_.space_ptr();

                const auto &sub_0 = space[0];

                const std::size_t n_quad_points = fe_[sub_0.subspace_id()]->get_phi()[0].size();

                unsigned int n_shape_functions = 0;
                space.each([&n_shape_functions, this](const int index, const LibMeshFunctionSpace &subspace) {
                    n_shape_functions += fe_[subspace.subspace_id()]->n_shape_functions();
                });

                ret.resize(n_shape_functions);
                for (std::size_t i = 0; i < n_shape_functions; ++i) {
                    ret[i].resize(n_quad_points);
                    for (auto &r : ret[i]) {
                        r.resize(space.n_subspaces());
                        r.zero();
                    }
                }

                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    int offset = 0;
                    space.each([&offset, this, &ret, qp](const int sub_index, const LibMeshFunctionSpace &s) {
                        const auto &fe = fe_[s.subspace_id()];
                        assert((static_cast<bool>(fe)));
                        const auto &fun = fe->get_phi();
                        uint n_shape_i = fe->n_shape_functions();

                        for (uint j = 0; j < n_shape_i; ++j) {
                            ret[offset++][qp](sub_index) = fun[j][qp];
                        }
                    });
                }
            }

            void grad_aux(JacobianType &ret) const {
                assert(!fe_.empty());

                auto &space = *fun_.space_ptr();

                const auto &sub_0 = space[0];
                const std::size_t n_quad_points = fe_[sub_0.subspace_id()]->get_phi()[0].size();
                const uint dim = space[0].mesh().spatial_dimension();

                uint n_shape_functions = 0;
                space.each([this, &n_shape_functions](const int, const LibMeshFunctionSpace &subspace) {
                    n_shape_functions += fe_[subspace.subspace_id()]->n_shape_functions();
                });

                ret.resize(n_shape_functions);
                for (std::size_t i = 0; i < n_shape_functions; ++i) {
                    ret[i].resize(n_quad_points);
                    // TensorValue is by default initialized to 0s
                }

                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    int offset = 0;
                    space.each(
                        [&offset, this, &ret, &space, qp, dim](const int sub_index, const LibMeshFunctionSpace &s) {
                            const auto &fe = fe_[s.subspace_id()];
                            const uint n_shape_i = fe->n_shape_functions();

                            for (uint j = 0; j < n_shape_i; ++j, offset++) {
                                const auto &grad = fe->get_dphi()[j][qp];

                                for (uint d = 0; d < dim; ++d) {
                                    ret[offset][qp](sub_index, d) = grad(d);
                                }
                            }
                        });
                }
            }

            void init_coeffs(const int element_idx) {
                auto space_ptr = fun_.space_ptr();
                const auto &sub_0 = space_ptr->subspace(0);
                const auto &mesh = sub_0.mesh();
                const auto &dof_map = sub_0.dof_map();

                const auto &elem_ptr = mesh.elem(element_idx);

                std::vector<libMesh::dof_id_type> prod_indices;
                std::vector<libMesh::dof_id_type> indices;

                space_ptr->each([&](const int sub_index, const LibMeshFunctionSpace &space) {
                    dof_map.dof_indices(elem_ptr, indices, space.subspace_id());
                    prod_indices.insert(prod_indices.end(), indices.begin(), indices.end());
                });

                const std::size_t n_indices = prod_indices.size();
                element_coeffs_.resize(n_indices);

                Read<UVector> r(coeff_);
                assert(coeff_.implementation().has_ghosts() || mpi_world_size() == 1);
                coeff_.get(prod_indices, element_coeffs_);
            }

            void init_function() {
                std::vector<std::vector<DenseVectorT>> g;
                fun_aux(g);

                const std::size_t n_shape_functions = g.size();
                const std::size_t n_quad_points = g[0].size();

                fun_values_.resize(n_quad_points);

                for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    fun_values_[qp].resize(fun_.codim());
                    fun_values_[qp].zero();

                    for (std::size_t i = 0; i < n_shape_functions; ++i) {
                        if (std::is_rvalue_reference<decltype(g)>::value) {
                            g[i][qp] *= element_coeffs_[i];
                            fun_values_[qp] += g[i][qp];
                        } else {
                            auto temp = g[i][qp];
                            temp *= element_coeffs_[i];
                            fun_values_[qp] += temp;
                        }
                    }
                }
            }

            void init_gradient() {}

            void init_divergence() {}
        };

        enum { Order = 1 };

        Interpolate(const UVector &coeff, Fun &fun) : data_(std::make_shared<Data>(coeff, fun)) {}

        inline const Fun &fun() const { return data_->fun_; }

        inline const UVector &coefficient() const { return data_->coeff_; }

        inline std::string get_class() const override { return "Interpolate"; }

        inline libMesh::FEBase &fe(const uint subspace_id) {
            assert(data_);
            assert(subspace_id < data_->fe_.size());
            assert(data_->fe_[subspace_id]);

            return *data_->fe_[subspace_id];
        }

    private:
        std::shared_ptr<Data> data_;
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_INTERPOLATE_HPP
