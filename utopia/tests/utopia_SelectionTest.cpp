#include "utopia.hpp"
#include "utopia_IsSubTree.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class SelectionTest {
    private:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using IndexSet = typename utopia::Traits<Vector>::IndexSet;
        using Comm = typename utopia::Traits<Vector>::Communicator;

        void vector_selection_test() {
            auto &&comm = Comm::get_default();
            const int n = comm.size() * 3;
            Vector v(layout(comm, 3, n), 0.);

            auto r = range(v);

            {
                Write<Vector> w_v(v);
                for (auto i = r.begin(); i < r.end(); ++i) {
                    v.set(i, i);
                }
            }

            IndexSet s;
            s.resize(2);
            {
                Write<IndexSet> w(s);
                // FIXME
                s[0] = r.begin();
                s[1] = r.end() % n;
            }

            Vector selection = select(v, s);
            auto s_r = range(selection);

            {
                Read<Vector> r_s(selection);
                utopia_test_assert(selection.get(s_r.begin()) == r.begin());
                utopia_test_assert(selection.get(s_r.begin() + 1) == (r.end() % n));
            }

            Scalar sum_v_s = sum(select(v, s));
            utopia_test_assert(sum_v_s >= r.begin());
        }

        void matrix_selection_test() {
            const SizeType n = mpi_world_size() * 3;
            Matrix m = sparse(n, n, n);
            auto rr = row_range(m);

            {
                Write<Matrix> w_m(m);
                for (auto i = rr.begin(); i < rr.end(); ++i) {
                    for (SizeType j = 0; j < n; ++j) {
                        m.set(i, j, i * n + j);
                    }
                }
            }

            IndexSet r_s;
            IndexSet c_s;

            r_s.resize(2);
            c_s.resize(2);

            {
                Write<IndexSet> wr(r_s), wc(c_s);
                // FIXME use set instead
                r_s[0] = rr.begin();
                r_s[1] = rr.begin() + 1;

                c_s[0] = 0;
                c_s[1] = 2;
            }

            Matrix selection = select(m, r_s, c_s);

            {
                auto s_r = row_range(selection);
                Read<Matrix> r_s(selection);
                utopia_test_assert(selection.get(s_r.begin(), 0) == rr.begin() * n);
                utopia_test_assert(selection.get(s_r.begin(), 1) == rr.begin() * n + 2);
            }

            Matrix row_selection = select(m, r_s, IndexSet() /* nvcc bug if we omit the default parameter */);

            {
                auto s_r = row_range(row_selection);
                Read<Matrix> r_s(row_selection);

                for (SizeType i = 0; i < n; ++i) {
                    utopia_test_assert(row_selection.get(s_r.begin(), i) == (rr.begin() * n + i));
                }
            }
        }

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

    public:
        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(vector_selection_test);
            UTOPIA_RUN_TEST(matrix_selection_test);
        }
    };

    static void selection(){
    // FIXME
#ifdef WITH_BLAS
    //     SelectionTest<BlasMatrixd, BlasVectord>().run();
#endif  // WITH_BLAS

#ifdef WITH_PETSC
    //   SelectionTest<PetscMatrix, PetscVector>().run();
#endif  // WITH_PETSC

        // #ifdef WITH_TRILINOS
        // FIME missing implementation in TpetraMatrix
        //         SelectionTest<TpetraMatrixd, TpetraVectord>().run();
        // #endif //WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(selection);
}  // namespace utopia
