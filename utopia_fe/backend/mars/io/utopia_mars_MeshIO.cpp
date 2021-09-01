// #include "utopia_mars_MeshIO.hpp"

// #include "utopia_mars_Mesh.hpp"

// namespace utopia {
//     namespace mars {

//         class MeshIO::Impl : public Configurable {
//         public:
//             void read(Input &in) override {}

//             bool write(const Path &write_path) {}

//             Impl(Mesh &mesh) : mesh(mesh), io_broker(utopia::make_unique<IOBroker>(mesh.comm().raw_comm())) {}

//         public:
//             Mesh &mesh;
//         };

//         void MeshIO::enable_interpolation_mode() { impl_->time_match = ::mars::io::MeshField::LINEAR_INTERPOLATION; }

//         void MeshIO::read(Input &in) { impl_->read(in); }
//         bool MeshIO::load() { return impl_->load(); }
//         bool MeshIO::write(const Path &write_path) { return impl_->write(write_path); }

//     }  // namespace mars
// }  // namespace utopia
