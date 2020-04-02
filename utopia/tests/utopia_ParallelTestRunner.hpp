#ifndef UTOPIA_PARALLEL_TEST_RUNNER_HPP
#define UTOPIA_PARALLEL_TEST_RUNNER_HPP

#include "utopia_Traits.hpp"
#include "utopia_Input.hpp"

namespace utopia {

    template<class Comm>
    class UnitTest {
    public:
        virtual ~UnitTest() {}
        virtual void set_up() {}
        virtual void tear_down() {}
        virtual void run() = 0;

        UnitTest() {}

        void set_comm(const Comm &comm)
        {
            comm_ = comm;
        }

        Comm &comm()
        {
            return comm_;
        }

        const Comm &comm() const
        {
            return comm_;
        }

    private:
        Comm comm_;
    };

    template<class Tensor>
    class AlgebraUnitTest : public UnitTest<typename Traits<Tensor>::Communicator>{
    public:
        virtual ~AlgebraUnitTest() {}
    };

    template<class Comm>
    class ParallelTestRunner {
    public:
        virtual ~ParallelTestRunner() {}

        void run(UnitTest<Comm> &test, const Comm &comm)
        {
            const int comm_size = comm.size();

            //1) we run on comm
            test.set_comm(comm);
            test.set_up();
            test.run();
            test.tear_down();

            //1) we create sub-comms
            for(auto &p : parititions_) {
                if(p < comm_size) {
                    Comm sub_comm = comm.split(comm.rank()/p);

                    test.set_comm(sub_comm);
                    test.set_up();
                    test.run();
                    test.tear_down();
                } else {
                    break;
                }
            }
        }

        ParallelTestRunner()
        : parititions_{{2, 4}}
        {}

    private:
        std::vector<int> parititions_;
    };

    template<class Comm>
    void run_parallel_test(
        UnitTest<Comm> &test,
        const Comm &comm = Comm()
    )
    {
        ParallelTestRunner<Comm>().run(test, comm);
    }

    template<class Test>
    void run_parallel_test()
    {
        Test test;
        run_parallel_test(test);
    }

}

#define UTOPIA_TEST_CLASS(macro_ClassName, macro_Matrix_, macro_Vector_) \
    template<class macro_Matrix_, class macro_Vector_> \
    class macro_ClassName final : public UnitTest< typename Traits<macro_Vector_>::Communicator >


#define UTOPIA_EXPOSE_TYPES(macro_Tensor_) \
    using Traits       = utopia::Traits<macro_Tensor_>; \
    using Scalar       = typename Traits::Scalar;   \
    using SizeType     = typename Traits::SizeType; \
    using IndexSet     = typename Traits::IndexSet; \
    using Comm         = typename Traits::Communicator; \
    using Layout       = typename Traits::Layout;   \
    using MatrixLayout = typename Traits::MatrixLayout;

#endif //UTOPIA_PARALLEL_TEST_RUNNER_HPP
