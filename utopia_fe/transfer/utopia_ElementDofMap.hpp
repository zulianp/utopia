#ifndef UTOPIA_ELEMENT_DOF_MAP_HPP
#define UTOPIA_ELEMENT_DOF_MAP_HPP

#include "moonolith_serializable.hpp"
#include "moonolith_input_stream.hpp"
#include "moonolith_output_stream.hpp"

#include <vector>

namespace utopia {

    template<typename T>
    inline void write_vector(
        const std::vector<T> &v,
        moonolith::OutputStream &os)
    {
        int n = v.size();
        os << n;
        os.write(&v[0], n);
    }

    template<typename T>
    inline void read_vector(
        std::vector<T> &v,
        moonolith::InputStream &is)
    {
        int n;
        is >> n;
        v.resize(n);
        is.read(&v[0], n);
    }

    class ElementDofMap : public moonolith::Serializable {
    public:
        ElementDofMap() : global_id(-1) {}

        inline void read(moonolith::InputStream &is) override
        {
            read_vector(global, is);
            is >> global_id;
        }

        inline void write(moonolith::OutputStream &os) const override
        {
            write_vector(global, os);
            os << global_id;
        }

        inline bool empty() const
        {
            return global.empty();
        }

        inline std::size_t n_dofs() const { return global.size(); }
        inline long dof(const std::size_t idx) const { assert(idx < global.size()); return global[idx]; }

        std::vector<long> global;
        long global_id;
    };
}

#endif //UTOPIA_ELEMENT_DOF_MAP_HPP
