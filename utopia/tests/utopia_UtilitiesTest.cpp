#include "utopia.hpp"
#include "utopia_AutoDiff.hpp"  //simplify_test
#include "utopia_Blocks.hpp"
// #include "utopia_CoreDecprecatedHeaders.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Eval_Blocks.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class BlockTest {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexSet = typename Traits::IndexSet;
        using Comm = typename Traits::Communicator;

        void run() { UTOPIA_RUN_TEST(block_test); }

    private:
        void block_test() {
            int n1 = 5;
            int n2 = 5;
            int n3 = 3;
            int n4 = 4;

            auto &&comm = Comm::get_default();

            auto vl1 = layout(comm, Traits::decide(), n1);
            auto vl2 = layout(comm, Traits::decide(), n2);
            auto vl3 = layout(comm, Traits::decide(), n3);
            auto vl4 = layout(comm, Traits::decide(), n4);

            auto id_ptr_11 = std::make_shared<Matrix>();
            id_ptr_11->identity(layout(comm, Traits::decide(), Traits::decide(), n1, n2), 1.0);
            auto id_ptr_12 = std::make_shared<Matrix>();
            id_ptr_12->identity(layout(comm, Traits::decide(), Traits::decide(), n3, n2), 2.0);
            auto id_ptr_22 = std::make_shared<Matrix>();
            id_ptr_22->identity(layout(comm, Traits::decide(), Traits::decide(), n3, n4), 3.0);

            Blocks<Matrix> b_mat(2, 2, {id_ptr_11, nullptr, id_ptr_12, id_ptr_22});

            auto s = size(b_mat);
            utopia_test_assert(s.get(0) == (n1 + n3));
            utopia_test_assert(s.get(1) == (n2 + n4));

            Matrix mat = b_mat;
            Vector twos(vl2, 2.);
            assert(twos.size() == n2);
            Vector ones(vl4, 1.);
            assert(ones.size() == n4);

            Vector vec = blocks(twos, ones);
            assert(vec.size() == (n2 + n4));
            Vector r = mat * vec;
            assert(r.size() == (n1 + n3));
            Vector r1(vl1, 0.0), r2(vl3, 0.0);
            assert(r1.size() == n1);
            assert(r2.size() == n3);

            undo_blocks(r, r1, r2);

            assert(r1.size() == n1);
            assert(r2.size() == n3);

            utopia_test_assert(approxeq(double(sum(r)), 2. * size(r1).get(0) + 7. * size(r2).get(0)));

            utopia_test_assert(approxeq(double(sum(r1)), 2. * size(r1).get(0)));

            utopia_test_assert(approxeq(double(sum(r2)), 7. * size(r2).get(0)));
        }
    };

    template <class Matrix, class Vector>
    class UtilitiesTest {
    private:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexSet = typename Traits::IndexSet;
        using Comm = typename Traits::Communicator;

        void csv_read_write() {
            Path path = Path(Utopia::instance().get("data_path")) / "csv/test.csv";
            CSV csv;

            utopia_test_assert(csv.read(path));
            utopia_test_assert(csv.write("./out.csv"));

            int val = -1.;
            csv.get(1, 0, val);
            utopia_test_assert(val == 0);
            csv.get(1, 1, val);
            utopia_test_assert(val == 1);
            csv.get(1, 2, val);
            utopia_test_assert(val == 2);

            csv.get(2, 0, val);
            utopia_test_assert(val == 1);
            csv.get(2, 1, val);
            utopia_test_assert(val == 2);
            csv.get(2, 2, val);
            utopia_test_assert(val == 3);
        }

        void factory_test() {
            Matrix m;
            m.identity(serial_layout(2, 2));
            auto size = m.size();
            utopia_test_assert(size.get(0) == 2);
            utopia_test_assert(size.get(1) == 2);

            m.read([](SizeType x, SizeType y, double entry) {
                if (x == y) {
                    utopia_test_assert(approxeq(1.0, entry));
                } else {
                    utopia_test_assert(approxeq(0.0, entry));
                }
            });

            Matrix m2;
            m2.dense(serial_layout(2, 2), -4);
            size = m2.size();
            utopia_test_assert(size.get(0) == 2);
            utopia_test_assert(size.get(1) == 2);

            m2.read([](SizeType, SizeType, double entry) { utopia_test_assert(approxeq(-4.0, entry)); });
        }

        void wrapper_test() {
            Vector v(serial_layout(2), 1.0);
            {
                Write<Vector> w(v);
                v.set(1, 2.0);
            }
            Matrix m;
            m.dense(serial_layout(2, 2), 0.0);
            {
                Write<Matrix> w(m);
                m.set(0, 0, 1.0);
                m.set(1, 1, 1.0);
            }

            Vector res;

            res = m * v * 0.1;

            {
                Read<Vector> w_res(res);
                utopia_test_assert(approxeq(0.1, res.get(0)));
                utopia_test_assert(approxeq(0.2, res.get(1)));
            }
        }

        // FIXME
        void range_test() {
            Matrix m1;
            m1.identity(serial_layout(3, 3));

            ////////////////////////////////////////////////////////////

            auto m1_view = view(m1, Range(0, 1), Range(0, 3));
            Matrix m2 = m1_view;

            m2.read([](SizeType /*x*/, SizeType y, double entry) {
                utopia_test_assert(approxeq(y == 0 ? 1.0 : 0.0, entry));
            });

#ifdef UTOPIA_WITH_PETSC
            // NOTE(eric): range assignment is NYI in Petsc backend
            if (std::is_same<Matrix, PetscMatrix>::value) {
                return;
            }
#endif

            ////////////////////////////////////////////////////////////

            Matrix m3 = m1;
            view(m3, Range(0, 1), Range(0, 3)) = m2;

            m3.read(
                [](SizeType x, SizeType y, double entry) { utopia_test_assert(approxeq(x == y ? 1.0 : 0.0, entry)); });

            Matrix diff = m1 - m3;
            diff.read([](SizeType /*x*/, SizeType /*y*/, double entry) { utopia_test_assert(approxeq(0.0, entry)); });
        }

        void factory_and_operations_test() {
            const int n = 3;
            Matrix m;
            m.dense(serial_layout(n, n), 0.0);

            {
                Write<Matrix> w(m);

                m.set(0, 0, 2.0);
                m.set(0, 1, -1.0);

                m.set(1, 0, -1.0);
                m.set(1, 1, 2.0);
                m.set(1, 2, -1.0);

                m.set(2, 1, -1.0);
                m.set(2, 2, -1.0);
            }
        }

        void variable_test() {
            using namespace std;

            auto v = make_shared<Vector>();
            v->values(serial_layout(20), 1.0);

            auto w = make_shared<Vector>();
            w->values(serial_layout(20), 3.0);

            Matrix m;
            m.identity(serial_layout(20, 20));

            Variable<Vector, 0> var_0(v);
            Variable<Vector, 1> var_1(w);

            auto expr = 0.5 * dot(m * var_0, var_1 - var_0);

            auto &var_0_found = find_variable<Vector, 0>(expr);
            utopia_test_assert(v.get() == &var_0_found.expr());

            auto &var_1_found = find_variable<Vector, 1>(expr);
            utopia_test_assert(w.get() == &var_1_found.expr());

            auto q = make_shared<Vector>();
            q->values(serial_layout(20), 4.0);

            var_0_found.set(w);
            var_1_found.set(q);

            double result = expr;
            utopia_test_assert(approxeq(10., result));
        }

    public:
        static void print_backend_info() {
            mpi_world_barrier();
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
            mpi_world_barrier();
        }

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(csv_read_write);
            UTOPIA_RUN_TEST(factory_test);
            UTOPIA_RUN_TEST(wrapper_test);
            UTOPIA_RUN_TEST(range_test);
            UTOPIA_RUN_TEST(factory_and_operations_test);
            UTOPIA_RUN_TEST(variable_test);
        }
    };

    static void authored_work_test() {
        std::stringstream ss;
        CitationsDB::instance().describe(ss);

        if (mpi_world_rank() == 0) {
            utopia_test_assert(!ss.str().empty());
        } else {
            utopia_test_assert(ss.str().empty());
        }
    }

    static void describe_test() {
        std::stringstream ss;

        Describable d;
        d.describe(ss);

        utopia_test_assert(!ss.str().empty());
    }

    static void describe_chrono_test() {
        std::stringstream ss;

        Chrono c;
        c.start();
        c.stop();
        c.describe(ss);

        utopia_test_assert(!ss.str().empty());

        c.rescale_duration(10);

        c.describe(ss);

        utopia_test_assert(!ss.str().empty());
    }

    static void utilities() {
        UTOPIA_RUN_TEST(authored_work_test);
        UTOPIA_RUN_TEST(describe_test);
        UTOPIA_RUN_TEST(describe_chrono_test);

#ifdef UTOPIA_WITH_BLAS
        UtilitiesTest<BlasMatrixd, BlasVectord>().run();
#endif  // UTOPIA_WITH_BLAS

#ifdef UTOPIA_WITH_PETSC
        BlockTest<PetscMatrix, PetscVector>().run();

        UtilitiesTest<PetscMatrix, PetscVector>().run();
        BlockTest<PetscMatrix, PetscVector>().run();

#endif  // UTOPIA_WITH_PETSC

        // FIXME
        if (mpi_world_size() == 1) {
#ifdef UTOPIA_WITH_TRILINOS
            BlockTest<TpetraMatrixd, TpetraVectord>().run();
#endif  // UTOPIA_WITH_TRILINOS
        }
    }

    UTOPIA_REGISTER_TEST_FUNCTION(utilities);
}  // namespace utopia
