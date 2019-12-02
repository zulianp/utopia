#include "utopia_LinePostProcessor.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_SideSetAssignment.hpp"
#include "MortarAssemble.hpp"

#include "libmesh/parallel_mesh.h"
#include "libmesh/elem.h"

#include <vector>

namespace utopia {

    template<class FunctionSpace, class Vector>
    class LinePostProcessor<FunctionSpace, Vector>::Impl : public PostProcessor<FunctionSpace, Vector> {
    public:
        void read(Input &in) override
        {
            in.get("p0-x", p0(0));
            in.get("p0-y", p0(1));
            in.get("p0-z", p0(2));

            in.get("p1-x", p1(0));
            in.get("p1-y", p1(1));
            in.get("p1-z", p1(2));

            in.get("n-samples", n_samples);

            generate_points();
        }

        void apply(FunctionSpace &V, const Vector &sol) override
        {
            auto &mesh = V.mesh();

            for(const auto &elem_ptr : mesh.active_local_element_ptr_range())
            {

            }
        }

        void apply(FunctionSpace &V, const Vector &sol, const Vector &other) override
        {

        }

        void describe(std::ostream &os) const override
        {

        }

        void export_values() const override
        {

        }

        Impl()
        : n_samples(10) {}

        class Sample {
        public:
            std::vector<double> values;
        };

    private:
        libMesh::Point p0, p1;
        std::size_t n_samples;

        std::vector<libMesh::Point> points;
        std::vector<std::string> sample_name;
        std::vector<Sample> sample;


        // moonolith::HPolytope<T, Dim> pol

        void generate_points()
        {
            double h = 1.0/(n_samples - 1);
            libMesh::Point v = h * (p1 - p0);

            points.resize(n_samples);
            points[0] = p0;

            for(std::size_t i = 1; i < n_samples; ++i)
            {
                points[i] = points[i-1] + v;
            }
        }
    };

    template<class FunctionSpace, class Vector>
    LinePostProcessor<FunctionSpace, Vector>::LinePostProcessor()
    : impl_(utopia::make_unique<Impl>())
    {}

    template<class FunctionSpace, class Vector>
    LinePostProcessor<FunctionSpace, Vector>::~LinePostProcessor()
    {}

    template<class FunctionSpace, class Vector>
    void LinePostProcessor<FunctionSpace, Vector>::apply(FunctionSpace &V, const Vector &sol)
    {
        impl_->apply(V, sol);
    }

    template<class FunctionSpace, class Vector>
    void LinePostProcessor<FunctionSpace, Vector>::apply(FunctionSpace &V, const Vector &sol, const Vector &other)
    {
        impl_->apply(V, sol, other);
    }

    template<class FunctionSpace, class Vector>
    void LinePostProcessor<FunctionSpace, Vector>::describe(std::ostream &os) const
    {
        impl_->describe(os);
    }

    template<class FunctionSpace, class Vector>
    void LinePostProcessor<FunctionSpace, Vector>::export_values() const
    {
        impl_->export_values();
    }

    template<class FunctionSpace, class Vector>
    void LinePostProcessor<FunctionSpace, Vector>::read(Input &in)
    {
        impl_->read(in);
    }

}
