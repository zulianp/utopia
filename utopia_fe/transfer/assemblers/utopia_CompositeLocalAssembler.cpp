#include "utopia_CompositeLocalAssembler.hpp"
#include "utopia.hpp"

namespace utopia {
    CompositeLocalAssembler::CompositeLocalAssembler(const std::vector<std::shared_ptr<LocalAssembler>> &assemblers)
    : assemblers_(assemblers)
    {}

    CompositeLocalAssembler::~CompositeLocalAssembler() {}

    bool CompositeLocalAssembler::assemble(
        const Elem &master,
        FEType master_type,
        const Elem &slave,
        FEType slave_type,
        Matrix &mat
        )
    {
        assert(false);
        return false;
    }

    bool CompositeLocalAssembler::assemble(
        const Elem &master,
        FEType master_type,
        const Elem &slave,
        FEType slave_type,
        std::vector<Matrix> &mat
        )
    {
        const std::size_t n = n_forms();
        mat.reserve(n);

        std::vector<Matrix> temp;
        for(std::size_t i = 0; i < n; ++i) {

            if(!assemblers_[i]->assemble(master, master_type, slave, slave_type, temp)) {
                return false;
            }

            for(auto &t : temp) {
                mat.push_back(t);
            }
        }

        return true;
    }

    int CompositeLocalAssembler::n_forms() const
    {
        int ret = 0;

        for(const auto &a_ptr : assemblers_) {
            ret += a_ptr->n_forms();
        }

        return ret;
    }

    CompositeLocalAssembler::Type CompositeLocalAssembler::type(const int index) const
    {
        int range_end = 0;
        int range_begin = 0;

        for(const auto &a_ptr : assemblers_) {
            range_begin = range_end;
            range_end += a_ptr->n_forms();

            if(index < range_end) {
                return a_ptr->type(index - range_begin);
            }
        }

        assert(false);
        m_utopia_error("should never arrive here");
        return MASTER_X_SLAVE;
    }

    void CompositeLocalAssembler::print_stats(std::ostream &os) const
    {
        for(const auto &a_ptr : assemblers_) {
            a_ptr->print_stats(os);
        }
    }
}
