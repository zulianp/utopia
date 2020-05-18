#include "utopia_SymbolicFunction.hpp"

#ifdef WITH_TINY_EXPR

#include "tinyexpr.h"
#include "utopia_make_unique.hpp"

#include <utility>
#include <vector>

namespace utopia {
    class SymbolicFunction::Impl {
    public:
        Impl(std::string expr)
            : x_(0.),
              y_(0.),
              z_(0.),
              t_(0.),
              expr_(std::move(expr)),
              vars_({
                  {"x", &x_, TE_VARIABLE, nullptr},
                  {"y", &y_, TE_VARIABLE, nullptr},
                  {"z", &z_, TE_VARIABLE, nullptr},
                  {"t", &t_, TE_VARIABLE, nullptr},
              }),
              err_(0) {
            /* Compile the expression with variables. */
            e_ = te_compile(expr_.c_str(), &vars_[0], vars_.size(), &err_);
        }

        inline bool valid() const { return err_ == 0; }

        ~Impl() { te_free(e_); }

        double eval() { return te_eval(e_); }

        inline const std::string &expr() const { return expr_; }

    public:
        double x_, y_, z_, t_;

    private:
        std::string expr_;

        std::vector<te_variable> vars_;
        int err_;

        te_expr *e_;
    };

    SymbolicFunction::~SymbolicFunction() = default;

    SymbolicFunction::SymbolicFunction(const std::string &expr) { impl_ = make_unique<Impl>(expr); }

    SymbolicFunction &SymbolicFunction::operator=(const SymbolicFunction &other) {
        if (this == &other) return *this;

        impl_ = make_unique<Impl>(other.impl_->expr());
        return *this;
    }

    SymbolicFunction::SymbolicFunction(const SymbolicFunction &other) {
        impl_ = make_unique<Impl>(other.impl_->expr());
    }

    bool SymbolicFunction::valid() const { return impl_->valid(); }

    double SymbolicFunction::eval(const double x) {
        impl_->x_ = x;
        impl_->y_ = 0.;
        impl_->z_ = 0.;
        impl_->t_ = 0.;
        return impl_->eval();
    }

    double SymbolicFunction::eval(const double x, const double y) {
        impl_->x_ = x;
        impl_->y_ = y;
        impl_->z_ = 0.;
        impl_->t_ = 0.;
        return impl_->eval();
    }

    double SymbolicFunction::eval(const double x, const double y, const double z) {
        impl_->x_ = x;
        impl_->y_ = y;
        impl_->z_ = z;
        impl_->t_ = 0.;
        return impl_->eval();
    }

    double SymbolicFunction::eval(const double x, const double y, const double z, const double t) {
        impl_->x_ = x;
        impl_->y_ = y;
        impl_->z_ = z;
        impl_->t_ = t;
        return impl_->eval();
    }

    double SymbolicFunction::eval(const std::vector<double> &x) {
        // FIXME
        assert(x.size() <= 4);

        switch (x.size()) {
            case 1: {
                return eval(x[0]);
            }

            case 2: {
                return eval(x[0], x[1]);
            }

            case 3: {
                return eval(x[0], x[1], x[2]);
            }

            case 4: {
                return eval(x[0], x[1], x[2], x[3]);
            }

            default: {
                std::cerr << "TODO" << std::endl;
                assert(false);
                return 0.;
            }
        }
    }
}  // namespace utopia

#endif  // WITH_TINY_EXPR
