#ifndef UTOPIA_UNIT_TEST_HPP
#define UTOPIA_UNIT_TEST_HPP

#include "utopia_TestAssert.hpp"

namespace utopia {

    template <class Comm>
    class UnitTest {
    public:
        virtual ~UnitTest() = default;
        virtual void set_up() {}
        virtual void tear_down() {}
        virtual void run() = 0;
        virtual void print_backend_info() const {}

        UnitTest() = default;

        void set_comm(const Comm &comm) { comm_ = comm; }

        Comm &comm() { return comm_; }

        const Comm &comm() const { return comm_; }

    private:
        Comm comm_;
    };

}  // namespace utopia

#endif  // UTOPIA_UNIT_TEST_HPP
