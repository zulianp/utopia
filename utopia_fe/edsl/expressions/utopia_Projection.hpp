#ifndef UTOPIA_FE_PROJECTION_HPP
#define UTOPIA_FE_PROJECTION_HPP

namespace utopia {
    template <class From, class To>
    class Projection : public Expression<Projection<From, To> > {
    public:
        Projection(const From &from, To &to) : from_(from), to_(to) {}

        inline const From &from() const { return from_; }

        inline To &to() const { return to_; }

    private:
        UTOPIA_STORE_CONST(From) from_;
        To &to_;
    };

    template <class From, class To>
    Projection<From, To> project(const From &from, To &to) {
        return Projection<From, To>(from, to);
    }
}  // namespace utopia

#endif  // UTOPIA_FE_PROJECTION_HPP
