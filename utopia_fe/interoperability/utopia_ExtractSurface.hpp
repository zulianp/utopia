#ifndef UTOPIA_EXTRACT_SURFACE_HPP
#define UTOPIA_EXTRACT_SURFACE_HPP

namespace utopia {
    template <class VolumeMesh, class SurfaceMesh>
    class ExtractSurface {};

    template <class VolumeMesh, class SurfaceMesh, typename... Args>
    inline void extract_surface(const VolumeMesh &in, SurfaceMesh &out, Args... args) {
        ExtractSurface<VolumeMesh, SurfaceMesh>::apply(in, out, args...);
    }
}  // namespace utopia

#endif  // UTOPIA_EXTRACT_SURFACE_HPP
