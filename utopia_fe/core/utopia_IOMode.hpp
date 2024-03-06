#ifndef UTOPIA_IO_MODE_HPP
#define UTOPIA_IO_MODE_HPP

namespace utopia {
    enum OutputMode {
        OUTPUT_MODE_OVERWRITE = 0,
        OUTPUT_MODE_RESTART = 1,
        OUTPUT_MODE_APPEND,
    };
}

#endif  // UTOPIA_IO_MODE_HPP