#include "utopia_DifferenceApp.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_libmesh.hpp"

#include "utopia_ui.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
// #include "libmesh/namebased_io.h"
#include "libmesh/distributed_mesh.h"
#include "utopia_MeshTransferOperator.hpp"

namespace utopia {

    class MeshData : public Configurable {
    public:
        MeshData(libMesh::Parallel::Communicator &comm) :
        mesh(std::make_shared<libMesh::DistributedMesh>(comm)),
        io(mesh.mesh()),
        space(make_ref(mesh))
        {}

        void read(Input &in) override
        {
            in.get("mesh", [this](Input &in) {
                std::string path;
                in.get("path", path);
                io.read(path.c_str());
            });

            mesh.mesh().prepare_for_use();

            in.get("space", space);

            init();
        }

        void init()
        {
            space.space().each([this](const int i, LibMeshFunctionSpace &s) {
                const auto &name = s.var_name();
                std::cout << "reading var \"" << name << "\"" << std::endl;
                io.copy_nodal_solution(
                                s.equation_system(),
                                name,
                                name,
                                1);
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

        Impl(libMesh::Parallel::Communicator &comm)
        : from_(comm), to_(comm), output_path_("diff.e")
        {}

        void read(Input &in) override
        {
            utopia::Utopia::instance().read(in);

            in.get("from", from_);
            in.get("to", to_);

            transfer_ = utopia::make_unique<MeshTransferOperator>(
                from_.mesh.mesh_ptr(),
                make_ref(from_.space.space()[0].dof_map()),
                to_.mesh.mesh_ptr(),
                make_ref(to_.space.space()[0].dof_map())
            );

            in.get("transfer", *transfer_);
            in.get("output-path", output_path_);
        }

        void run()
        {
            if(!transfer_->assemble()) {
                assert(false);
                return;
            }

            transfer_->apply(from_.data, diff_);

            diff_ -= to_.data;
            diff_ = abs(diff_);

            /////////////////////////////////////////////////////////////////

            std::cout << "norm1      : " << double(norm1(diff_))      << std::endl;
            std::cout << "norm2      : " << double(norm2(diff_))      << std::endl;
            std::cout << "norm_infty : " << double(norm_infty(diff_)) << std::endl;

            auto &V = to_.space.space();

            /////////////////////////////////////////////////////////////////

            double l2_norm = 0.0;

            auto u = trial(V);
            auto x = interpolate(diff_, u);
            utopia::assemble(inner(x, x) * dX, l2_norm);
            l2_norm = std::sqrt(l2_norm);

            std::cout << "l2_norm    : " << l2_norm << std::endl;

            /////////////////////////////////////////////////////////////////

            write(output_path_, to_.space.space()[0], diff_);
        }

    private:
        MeshData from_, to_;
        std::unique_ptr<MeshTransferOperator> transfer_;
        std::string output_path_;
        UVector diff_;
    };

    DifferenceApp::DifferenceApp()
    {}

    DifferenceApp::~DifferenceApp()
    {}

    void DifferenceApp::run(Input &in)
    {
        auto impl = utopia::make_unique<Impl>(this->comm());
        impl->read(in);
        impl->run();
    }

}