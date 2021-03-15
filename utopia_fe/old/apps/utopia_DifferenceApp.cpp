#include "utopia_DifferenceApp.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_ui.hpp"
// #include "libmesh/namebased_io.h"
#include "libmesh/distributed_mesh.h"
#include "utopia_Deformation.hpp"
#include "utopia_MeshTransferOperator.hpp"

namespace utopia {

    class MeshData : public Configurable {
    public:
        MeshData(libMesh::Parallel::Communicator &comm)
            : mesh(std::make_shared<libMesh::DistributedMesh>(comm)), io(mesh.mesh()), space(make_ref(mesh)) {}

        void read(Input &in) override {
            in.get("mesh", [this](Input &in) {
                std::string path;
                in.get("path", path);
                io.read(path.c_str());
            });

            mesh.mesh().prepare_for_use();

            in.get("space", space);

            init();
        }

        void init() {
            space.space().each([this](const int i, LibMeshFunctionSpace &s) {
                const auto &name = s.var_name();
                std::cout << "reading var \"" << name << "\"" << std::endl;
                io.copy_nodal_solution(s.equation_system(), name, name, 1);
            });

            utopia::convert(*space.space()[0].equation_system().solution, data);
        }

        UIMesh<libMesh::DistributedMesh> mesh;
        libMesh::ExodusII_IO io;
        UIFunctionSpace<LibMeshFunctionSpace> space;
        UVector data;
    };

    class DifferenceApp::Impl : public Configurable {
    public:
        class Diff {
        public:
            std::size_t n_elements;
            std::size_t n_nodes;
            double relative_l2_squared_norm;
            double l2_squared_norm;
            double ref_min;
            double ref_max;
            double delta;
            double vol;
            double norm_1;
            double norm_2;
            double norm_infty;

            bool write(const Path &path) const {
                std::ofstream os(path.c_str());

                if (os.good()) {
                    describe(os);
                    os.close();
                    return true;
                    ;
                }

                os.close();
                return false;
            }

            void describe(std::ostream &os) const {
                os << "\"n_elements\",\"n_nodes\",\"relative_l2_squared_norm\",\"l2_squared_norm\",\"ref_min\",\"ref_"
                      "max\",\"delta\",\"vol\",";
                os << "\"norm_1\",\"norm_2\",\"norm_infty\"\n";

                os << n_elements << "," << n_nodes << "," << relative_l2_squared_norm << "," << l2_squared_norm << ","
                   << ref_min << "," << ref_max << "," << delta << "," << vol << "," << norm_1 << "," << norm_2 << ","
                   << norm_infty << "\n";
            }

            void describe_pretty(std::ostream &os) const {
                os << "n_elements               : " << n_elements << std::endl;
                os << "n_nodes                  : " << n_nodes << std::endl;
                os << "relative_l2_squared_norm : " << relative_l2_squared_norm << std::endl;
                os << "l2_squared_norm          : " << l2_squared_norm << std::endl;
                os << "ref_min                  : " << ref_min << std::endl;
                os << "ref_max                  : " << ref_max << std::endl;
                os << "delta                    : " << delta << std::endl;
                os << "vol                      : " << vol << std::endl;
                os << "norm_1                   : " << norm_1 << std::endl;
                os << "norm_2                   : " << norm_2 << std::endl;
                os << "norm_infty               : " << norm_infty << std::endl;
            }
        };

        Impl(libMesh::Parallel::Communicator &comm)
            : from_(comm), to_(comm), output_path_("diff.e"), csv_ouput_path_("diff.csv") {}

        void read(Input &in) override {
            utopia::Utopia::instance().read(in);

            in.get("from", from_);
            in.get("to", to_);

            transfer_ = utopia::make_unique<MeshTransferOperator>(from_.mesh.mesh_ptr(),
                                                                  make_ref(from_.space.space()[0].dof_map()),
                                                                  to_.mesh.mesh_ptr(),
                                                                  make_ref(to_.space.space()[0].dof_map()));

            in.get("transfer", *transfer_);
            in.get("output-path", output_path_);
            in.get("csv-ouput-path", csv_ouput_path_);
        }

        void run() {
            if (!transfer_->assemble()) {
                assert(false);
                return;
            }

            transfer_->apply(from_.data, diff_);

            diff_ -= to_.data;
            diff_ = abs(diff_);

            /////////////////////////////////////////////////////////////////

            auto &V = to_.space.space();

            /////////////////////////////////////////////////////////////////

            double l2_squared_norm = 0.0;

            auto u = trial(V);
            auto x = interpolate(diff_, u);
            utopia::assemble(inner(x, x) * dX, l2_squared_norm);

            double ref_min = min(to_.data);
            double ref_max = max(to_.data);
            double delta = ref_max - ref_min;

            USparseMatrix mass_matrix;  // FIXME
            utopia::assemble(inner(trial(V), test(V)) * dX, mass_matrix);
            double vol = sum(mass_matrix);

            Diff diff;
            diff.n_elements = from_.mesh.mesh().n_active_elem();
            diff.n_nodes = from_.mesh.mesh().n_nodes();
            diff.relative_l2_squared_norm = l2_squared_norm / (vol * delta * delta);
            diff.l2_squared_norm = l2_squared_norm;
            diff.ref_min = ref_min;
            diff.ref_max = ref_max;
            diff.delta = delta;
            diff.vol = vol;
            diff.norm_1 = norm1(diff_);
            diff.norm_2 = norm2(diff_);
            diff.norm_infty = norm_infty(diff_);

            if (mass_matrix.comm().rank() == 0) {
                diff.describe_pretty(std::cout);
                diff.write(csv_ouput_path_);
            }

            /////////////////////////////////////////////////////////////////

            write(output_path_, to_.space.space()[0], diff_);
        }

    private:
        MeshData from_, to_;
        std::unique_ptr<MeshTransferOperator> transfer_;
        std::string output_path_, csv_ouput_path_;
        UVector diff_;
    };

    DifferenceApp::DifferenceApp() {}

    DifferenceApp::~DifferenceApp() {}

    void DifferenceApp::run(Input &in) {
        auto impl = utopia::make_unique<Impl>(this->comm());
        impl->read(in);
        impl->run();
    }

}  // namespace utopia