#ifndef UTOPIA_FE_TRANSFER_APP_HPP
#define UTOPIA_FE_TRANSFER_APP_HPP

#include "utopia_fe_base.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_FEModelFunction.hpp"
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
    class FETransferTimeSeriesApp : public Configurable {
    public:
        // Extract front-end associated objects
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using IO_t = utopia::IO<FunctionSpace>;

        void read(Input &in) override {
            in.get("verbose", verbose);

            if (verbose) {
                utopia::out() << "Reading from!\n";
            }

            input = utopia::make_unique<IO_t>(from_space);
            input->import_all_field_data(true);
            input->enable_interpolation_mode();

            in.get("from", [this](Input &in) {
                bool valid = input->open_input(in);
                assert(valid);
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
                    IO<FunctionSpace> input_db(to_space);
                    input_db.set_import_all_data(true);
                    input_db.open_input(in);

                    double time = 0;
                    in.get("load_time", time);
                    input_db.load_time_step(time);

                    Field<FunctionSpace> displacement_field{"disp"};
                    displacement_field.set_space(make_ref(to_space));

                    int sd = to_space.mesh().spatial_dimension();
                    displacement_field.set_tensor_size(sd);

                    auto vlo =
                        layout(to_space.comm(), to_space.mesh().n_local_nodes() * sd, to_space.mesh().n_nodes() * sd);
                    displacement_field.set_data(std::make_shared<Vector_t>(vlo, 0));
                    input_db.read_nodal(displacement_field, true);

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
            in.get("export_example_coupled_system", export_example_coupled_system);

            // if (export_operator_imbalance && from_space.comm().size() != 1) {
            //     if (!from_space.comm().rank()) {
            //         utopia::err() << "Option \"export_operator_imbalance : true\" only works for serial runs!\n";
            //     }

            //     // Utopia::Abort();
            // }

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

            if (!transfer.init(make_ref(from_space), make_ref(to_space))) {
                return;
            }

            if (verbose) {
                from_space.comm().root_print("Exiting transfer!\n", utopia::out().stream());
            }

            Field<FunctionSpace> to_field;
            to_space.create_field(to_field);

            IO<FunctionSpace> output(to_space);
            output.set_output_path(output_path);

            int step = 0;
            for (auto t : input->get_time_steps()) {
                utopia::out() << "Resampling time step t = " << t << "\n";
                input->load_time_step(t);
                input->read_nodal(field);
                transfer.apply(field.data(), to_field.data());
                output.write(to_field.data(), step++, t);
            }
        }

        FETransferTimeSeriesApp() {}

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
        bool export_example_coupled_system{false};
        Scalar_t rescale_imbalance{1};

        std::shared_ptr<IO_t> input;
    };

}  // namespace utopia

#endif  // UTOPIA_FE_TRANSFER_APP_HPP
