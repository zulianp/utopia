#ifndef UTOPIA_SIDE_SET_ASSIGNMENT_HPP
#define UTOPIA_SIDE_SET_ASSIGNMENT_HPP

#include <memory>
#include "utopia_Input.hpp"

namespace utopia {

    template <class Mesh>
    class SideSetAssignment : public Configurable {
    public:
        void read(Input &in) override;
        void apply(Mesh &mesh);

        SideSetAssignment();
        ~SideSetAssignment();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };
}  // namespace utopia

#endif  // UTOPIA_SIDE_SET_ASSIGNMENT_HPP
