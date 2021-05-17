#ifndef UTOPIA_MESH_TEST_HPP
#define UTOPIA_MESH_TEST_HPP

#include "utopia_Traits.hpp"

#include "utopia_Testing.hpp"
#include "utopia_UnitTest.hpp"

namespace utopia {

    template <class Mesh>
    class MeshTest final : public UnitTest<typename Traits<Mesh>::Communicator> {
    public:
        using SizeType = typename Traits<Mesh>::SizeType;
        void unit_cube() {
            const SizeType n = this->comm().size() * 2;

            Mesh mesh(this->comm());
            mesh.unit_cube(n, n, n);

            const SizeType n_nodes = mesh.n_nodes();
            UTOPIA_TEST_EQ(n_nodes, ((n + 1) * (n + 1) * (n + 1)));
        }

        void read_write() {
            const SizeType n = this->comm().size() * 2;

            Mesh mesh(this->comm());
            mesh.unit_cube(n, n, n);
            UTOPIA_TEST_TRUE(mesh.write("./mesh_test_read_write.e"));
            UTOPIA_TEST_TRUE(mesh.read("./mesh_test_read_write.e"));
        }

        void run() override {
            UTOPIA_RUN_TEST(unit_cube);
            UTOPIA_RUN_TEST(read_write);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MESH_TEST_HPP