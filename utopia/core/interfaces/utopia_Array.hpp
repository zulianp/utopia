#include "utopia_Enums.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Range.hpp"

namespace utopia {

    template <typename EntryType, typename SizeType_>
    class Array {
    public:
        using EntryType = EntryType;
        using SizeType = SizeType_;

        virtual ~Array() {}

        // locks
        virtual void read_lock() = 0;
        virtual void write_lock(WriteMode) = 0;

        virtual void read_unlock() = 0;
        virtual void write_unlock(WriteMode) = 0;

        // basic mutators
        virtual void set(const SizeType &i, const EntryType &value) = 0;
        virtual EntryType get(const SizeType &i) = 0;

        // print function
        virtual void describe() const = 0;

        // utility functions
        virtual bool empty() const = 0;
        virtual void clear() = 0;
        virtual SizeType size() const = 0;
    };

}  // namespace utopia
