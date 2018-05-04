/*
* @Author: Eric Botter
* @Date:   2016-11-07
*/
#include "utopia.hpp"
#include "utopia_SelectionTest.hpp"
// #include "utopia_Collection.hpp"
#include "utopia_IsSubTree.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class SelectionTest {
    private:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;

        void vector_selection_test()
        {
            typedef typename utopia::Traits<Vector>::SizeType SizeType;

            const int n = mpi_world_size() * 3;
            Vector v = zeros(n);
            auto r = range(v);

            {
                Write<Vector> w_v(v);
                for(auto i = r.begin(); i < r.end(); ++i) {
                    v.set(i, i);
                }
            }

            std::vector<SizeType> s;
            s.push_back(r.begin());
            s.push_back(r.end() % n);

            Vector selection = v.select(s);
            auto s_r = range(selection);

            {
                Read<Vector> r_s(selection);
                utopia_test_assert(selection.get(s_r.begin()) == r.begin());
                utopia_test_assert(selection.get(s_r.begin() + 1) == (r.end() % n));
            }

            Scalar sum_v_s = sum(v.select(s));
            utopia_test_assert(sum_v_s >= r.begin());
        }

        void matrix_selection_test()
        {
            typedef typename utopia::Traits<Vector>::SizeType SizeType;
            
            const SizeType n = mpi_world_size() * 3;
            Matrix m = zeros(n, n);
            auto rr = row_range(m);

            {
                Write<Matrix> w_m(m);
                for(auto i = rr.begin(); i < rr.end(); ++i) {
                    for(SizeType j = 0; j < n; ++j) {
                        m.set(i, j, i * n + j);
                    }
                }
            }

            std::vector<SizeType> r_s;
            std::vector<SizeType> c_s;

            r_s.push_back(rr.begin());
            r_s.push_back(rr.begin() + 1);

            c_s.push_back(0);
            c_s.push_back(2);

            Matrix selection = m.select(r_s, c_s);

            {
                auto s_r = row_range(selection);
                Read<Matrix> r_s(selection);
                utopia_test_assert(selection.get(s_r.begin(), 0) == rr.begin() * n);
                utopia_test_assert(selection.get(s_r.begin(), 1) == rr.begin() * n + 2);
            }

            Matrix row_selection = m.select(r_s);

            {
                auto s_r = row_range(row_selection);
                Read<Matrix> r_s(row_selection);

                for(SizeType i = 0; i < n; ++i) {
                    utopia_test_assert(row_selection.get(s_r.begin(), i) == (rr.begin() * n + i));
                }
            }

        }

        static void print_backend_info()
        {
            if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

    public:
        void run()
        {
            print_backend_info();
            UTOPIA_RUN_TEST(vector_selection_test);
            UTOPIA_RUN_TEST(matrix_selection_test);
        }
    };

    void run_selection_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("SelectionTest");

#ifdef WITH_BLAS
        SelectionTest<Matrixd, Vectord>().run();
#endif //WITH_BLAS

#ifdef WITH_PETSC
        SelectionTest<DMatrixd, DVectord>().run();
#endif //WITH_PETSC

// #ifdef WITH_TRILINOS
//         SelectionTest<TSMatrixd, TVectord>().run();
// #endif //WITH_TRILINOS

        UTOPIA_UNIT_TEST_END("SelectionTest");
    }
}
