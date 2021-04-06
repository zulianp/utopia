#ifndef UTOPIA_STK_IO_HPP
#define UTOPIA_STK_IO_HPP

#include "utopia_Input.hpp"
#include "utopia_Path.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"

#include <memory>
#include <string>

namespace utopia {
    namespace stk {

        class IO final : public Configurable {
        public:
            void read(Input &in) override;
            bool load();
            bool write(const Path &write_path);

            void set_read_path(const Path &path);
            void set_read_specification(const std::string &format);

            IO(Mesh &mesh);
            ~IO();

        public:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_IO_HPP
