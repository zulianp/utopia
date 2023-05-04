#ifndef UTOPIA_FE_TRANSFER_APP_HPP
#define UTOPIA_FE_TRANSFER_APP_HPP

#include "utopia_fe_base.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_Field.hpp"
#include "utopia_ILU.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

#include <memory>

namespace utopia {

    template <class FunctionSpace>
    class FETransferApp : public Configurable {
    public:
        // Extract front-end associated objects
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;

        void read(Input &in) override {
            in.get("verbose", verbose);

            if (verbose) {
                utopia::out() << "Reading from!\n";
            }

            in.get("from", [this](Input &in) {
                bool read_state = false;
                in.get("read_state", read_state);

                if (read_state) {
                    from_space.read_with_state(in, field);

                    const Scalar_t norm_field = norm2(field.data());
                    const Size_t size_field = field.data().size();

                    const auto n_elements = from_space.mesh().n_elements();
                    const auto n_nodes = from_space.mesh().n_nodes();

                    if (from_space.comm().rank() == 0) {
                        utopia::out() << "[From Mesh] n elements: " << n_elements << ", n nodes: " << n_nodes << '\n';

                        utopia::out() << "[From Field] size: " << size_field << ", norm: " << norm_field << '\n';
                    }

                } else {
                    from_space.read(in);
                }
            });

            if (from_space.empty()) {
                return;
            }

            if (verbose) {
                from_space.comm().root_print("Reading to!\n", utopia::out().stream());
            }

            // in.get("to", to_space);
            in.get("to", [this](Input &in) {
                bool displace = false;
                in.get("displace", displace);

                if (displace) {
                    Field<FunctionSpace> displacement_field;
                    to_space.read_with_state(in, displacement_field);
                    to_space.displace(displacement_field.data());

                    const Scalar_t norm_displacement_field = norm2(displacement_field.data());
                    const Size_t size_displacement_field = displacement_field.data().size();

                    const auto n_elements = to_space.mesh().n_elements();
                    const auto n_nodes = to_space.mesh().n_nodes();

                    if (to_space.comm().rank() == 0) {
                        utopia::out() << "[To Mesh] n elements: " << n_elements << ", n nodes: " << n_nodes << '\n';
                        utopia::out() << "[To Field] size: " << size_displacement_field
                                      << ", norm: " << norm_displacement_field << '\n';
                    }
                } else {
                    to_space.read(in);
                }
            });

            if (to_space.empty()) {
                return;
            }

            if (field.empty()) {
                from_space.create_field(field);
                field.data().set(1.0);
            }

            in.get("transfer", transfer);
            in.get("output_path", output_path);
            in.get("export_from_function", export_from_function);
            in.get("export_operator_imbalance", export_operator_imbalance);
            in.get("rescale_imbalance", rescale_imbalance);

            if (export_operator_imbalance && from_space.comm().size() != 1) {
                if (!from_space.comm().rank()) {
                    utopia::err() << "Option \"export_operator_imbalance : true\" only works for serial runs!\n";
                }

                Utopia::Abort();
            }

            if (verbose) {
                from_space.comm().root_print("Exiting read!\n", utopia::out().stream());
            }
        }

        void run() {
            if (export_from_function) {
                from_space.write("./from_out.e", field.data());

                if (verbose) {
                    from_space.comm().root_print("Exported from function!\n", utopia::out().stream());
                }
            }

            // Matrix_t n2e;
            // from_space.create_node_to_element_matrix(n2e);

            if (!transfer.init(make_ref(from_space), make_ref(to_space))) {
                return;
            }

            if (verbose) {
                from_space.comm().root_print("Exiting transfer!\n", utopia::out().stream());
            }

            Field<FunctionSpace> to_field;
            to_space.create_field(to_field);

            transfer.apply(field.data(), to_field.data());
            to_space.write(output_path, to_field.data());

            if (export_operator_imbalance) {
                auto mat = transfer.transfer_matrix();

                mat->transform(
                    UTOPIA_LAMBDA(const Size_t &, const Size_t &, const Scalar_t &v)->Scalar_t { return 1; });

                // Compute imbalance = T^T * T * 1
                Vector_t ones(col_layout(*mat), 1);
                Vector_t T_ones = (*mat) * ones;
                Vector_t imbalance = transpose(*mat) * T_ones;

                from_space.write("imbalance.e", imbalance);

                Matrix_t n2e;
                from_space.create_node_to_element_matrix(n2e);

                Vector_t en_vec = n2e * imbalance;

                auto enl = layout(en_vec);

                int nnodexelement = enl.size() / from_space.mesh().n_elements();

                Vector_t e_vec(
                    layout(enl.comm(), from_space.mesh().n_local_elements(), from_space.mesh().n_elements()));

                {
                    auto en_view = const_local_view_device(en_vec);
                    auto e_view = local_view_device(e_vec);

                    parallel_for(
                        local_range_device(e_vec), UTOPIA_LAMBDA(const Size_t i) {
                            //
                            Scalar_t val = 0;
                            for (int d = 0; d < nnodexelement; d++) {
                                val += en_view.get(i * nnodexelement + d);
                            }

                            e_view.set(i, val);
                        });
                }

                Scalar_t max_imbalance = max(e_vec);
                max_imbalance = std::max(max_imbalance, Scalar_t(1));
                // Normalize
                e_vec *= rescale_imbalance / max_imbalance;
                e_vec.shift(1.0);

                Field<FunctionSpace> elemental_field("cost", make_ref(from_space), make_ref(e_vec));
                from_space.backend_set_elemental_field(elemental_field);

                IO<FunctionSpace> io(from_space);
                io.set_output_path("cost.e");
                io.register_output_field("cost");
                io.write(1, 1);
            }
        }

        FETransferApp() {}

        bool is_valid() const { return !from_space.empty() && !to_space.empty() && !field.empty(); }

    private:
        FunctionSpace from_space;
        FunctionSpace to_space;
        Field<FunctionSpace> field;
        FETransfer<FunctionSpace> transfer;
        Path output_path{"./out.e"};
        bool export_from_function{false};
        bool verbose{true};
        bool export_operator_imbalance{false};
        Scalar_t rescale_imbalance{1};
    };

}  // namespace utopia

#endif  // UTOPIA_FE_TRANSFER_APP_HPP
