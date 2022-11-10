
#include "Box.hpp"
#include <assert.h>

namespace utopia {

    void get_row(const int row,
                 const libMesh::DenseMatrix<libMesh::Real> &mat,
                 libMesh::DenseVector<libMesh::Real> &result) {
        result.resize(mat.n());
        for (int i = 0; i < mat.n(); ++i) {
            result(i) = mat(row, i);
        }
    }

    Box::Box(const int n) : min_(n), max_(n) { reset(); }

    Box::Box() {}

    Box::~Box() {}

    void Box::reset(const int n) {
        min_.resize(n);
        max_.resize(n);
        reset();
    }

    void Box::reset() {
        std::fill(min_.get_values().begin(), min_.get_values().end(), std::numeric_limits<libMesh::Real>::max());
        std::fill(max_.get_values().begin(), max_.get_values().end(), -std::numeric_limits<libMesh::Real>::max());
    }

    Box &Box::operator+=(const libMesh::Point &point) {
        using std::max;
        using std::min;

        int n = get_dims();
        assert(n > 0);

        for (int i = 0; i < n; ++i) {
            min_(i) = min(point(i), min_(i));
            max_(i) = max(point(i), max_(i));
        }

        return *this;
    }

    Box &Box::operator+=(const Box &box) {
        using std::max;
        using std::min;

        if (empty()) {
            *this = box;
            return *this;
        }

        int n = get_dims();

        for (int i = 0; i < n; ++i) {
            min_(i) = min(min_(i), box.min_(i));
            max_(i) = max(max_(i), box.max_(i));
        }

        return *this;
    }

    bool Box::intersects(const Box &other) const {
        int n = get_dims();

        assert(n == other.get_dims() && "must have same dimensions");

        for (int i = 0; i < n; ++i) {
            if (other.get_max(i) < get_min(i) || get_max(i) < other.get_min(i)) {
                return false;
            }
        }

        return true;
    }

    bool Box::intersects(const Box &other, const libMesh::Real tol) const {
        int n = get_dims();

        assert(n == other.get_dims() && "must have same dimensions");

        for (int i = 0; i < n; ++i) {
            if (other.get_max(i) + tol < get_min(i) || get_max(i) + tol < other.get_min(i)) {
                return false;
            }
        }

        return true;
    }

    void Box::enlarge(const libMesh::Real value) {
        auto abs_value = std::abs(value);
        const int n = get_dims();
        for (int i = 0; i < n; ++i) {
            min_(i) -= abs_value;
            max_(i) += abs_value;
        }
    }

    void Box::print(std::ostream &os) const {
        os << "[\n";
        for (int i = 0; i < get_dims(); ++i) {
            os << "\t" << get_min(i) << ", " << get_max(i) << "\n";
        }
        os << "]\n";
    }

}  // namespace utopia
