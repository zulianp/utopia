#ifndef UTOPIA_INTREPID2_ARBORX_INTERSECTOR_HPP
#define UTOPIA_INTREPID2_ARBORX_INTERSECTOR_HPP

#include <memory>

namespace utopia {

    template <typename CellPoints>
    class Intrepid2ArborXIntersector {
    public:
        Intrepid2ArborXIntersector();
        ~Intrepid2ArborXIntersector();

        bool detect(const CellPoints &points1, const CellPoints &points2);

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_ARBORX_INTERSECTOR_HPP
