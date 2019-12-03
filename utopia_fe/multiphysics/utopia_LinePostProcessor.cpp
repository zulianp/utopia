#include "utopia_LinePostProcessor.hpp"
#include "utopia_make_unique.hpp"
#include "moonolith_interpolator_impl.hpp"
#include "moonolith_elem_node1.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_SideSetAssignment.hpp"
#include "MortarAssemble.hpp"

#include "utopia_LibMeshToMoonolithConvertions.hpp"

#include "libmesh/parallel_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"

#include <vector>

namespace utopia {

    template<int Dim>
    void make_element(const libMesh::Point &in, moonolith::Node1<double, Dim> &out)
    {
        for(int d = 0; d < Dim; ++d) {
            out.node(0)[d] = in(d);
        }
    }

    template<int Rows, int Cols>
    void convert(const moonolith::Matrix<double, Rows, Cols> &in, USerialMatrix &out)
    {
        const std::size_t rows = in.rows();
        const std::size_t cols = in.cols();
        out.resize(rows, cols);

        for(std::size_t i = 0; i < rows; ++i) {
            for(std::size_t j = 0; j < cols; ++j) {
                out.set(i, j, in(i, j));
            }
        }
    }

    template<class FunctionSpace, class Vector>
    class LinePostProcessor<FunctionSpace, Vector>::Impl : public PostProcessor<FunctionSpace, Vector> {
    public:
        void read(Input &in) override
        {
            std::cout << "LinePostProcessor::read(...)" << std::endl;

            in.get("p0-x", p0(0));
            in.get("p0-y", p0(1));
            in.get("p0-z", p0(2));

            in.get("p1-x", p1(0));
            in.get("p1-y", p1(1));
            in.get("p1-z", p1(2));

            int temp = n_samples;
            in.get("n-samples", temp);
            n_samples = temp;

            in.get("output-path", output_path);

            std::cout << "p0-x       :" << p0(0) << std::endl;
            std::cout << "p0-y       :" << p0(1) << std::endl;
            std::cout << "p0-z       :" << p0(2) << std::endl;
            std::cout << "p1-x       :" << p1(0) << std::endl;
            std::cout << "p1-y       :" << p1(1) << std::endl;
            std::cout << "p1-z       :" << p1(2) << std::endl;
            std::cout << "n-samples  :" << n_samples << std::endl;
            std::cout << "output-path:" << output_path << std::endl;
            generate_points();
        }

        void apply(FunctionSpace &V, const Vector &sol) override
        {
            using namespace moonolith;

            auto &mesh = V.mesh();

            UIndexSet ghost_nodes;
            convert(V.dof_map().get_send_list(), ghost_nodes);
            Vector g_sol = ghosted(V.dof_map().n_local_dofs(), V.dof_map().n_dofs(), ghost_nodes);
            g_sol = sol;
            synchronize(g_sol);

            const int dim = mesh.mesh_dimension();
            assert(dim == mesh.spatial_dimension());

            Node1<double, 2> node2;
            Node1<double, 3> node3;

            std::shared_ptr< Elem<double, 2, 2> > elem2;
            std::shared_ptr< Elem<double, 3, 3> > elem3;

            Interpolator<Elem<double, 2, 2>, Elem<double, 0, 2>> interp2;
            Interpolator<Elem<double, 3, 3>, Elem<double, 0, 3>> interp3;

            libMesh::FEType fe_type = V.type();

            std::vector<libMesh::dof_id_type> indices;
            USerialVector values, interp_value;
            USerialMatrix weights;

            sample.resize(n_samples);
            sample_name.resize(1);

            sample_name[0] = V.var_name();

            //FIXME use AABB to make it faster

            UIndexSet u_indices;
            for(const auto &elem_ptr : mesh.active_local_element_ptr_range())
            {

                V.dofs(*elem_ptr, indices);
                convert(indices, u_indices);

                g_sol.get(u_indices, values.entries());

                if(dim == 2) {
                    make<2>(*elem_ptr, fe_type, elem2);
                } else {
                    make<3>(*elem_ptr, fe_type, elem3);
                }

                for(std::size_t i = 0; i < n_samples; ++i) {

                    bool found = false;
                    if(dim == 2) {
                        make_element(points[i], node2);

                        if(interp2.assemble(*elem2, node2)) {
                            convert(interp2.coupling_matrix(), weights);
                            found = true;
                        }

                    } else {
                        assert(dim == 3);
                        make_element(points[i], node3);

                        if(interp3.assemble(*elem3, node3)) {
                            convert(interp3.coupling_matrix(), weights);
                            found = true;
                        }
                    }

                    if(found) {
                        interp_value = weights * values;
                        sample[i].values[0] = interp_value.get(0);
                    }
                }
            }

            export_values();
        }

        void apply(FunctionSpace &V, const Vector &sol, const Vector &other) override
        {
            this->apply(V, sol);
        }

        void describe(std::ostream &os) const override
        {
            os << "\"arc_length\",\"Points:0\",\"Points:1\",\"Points:2\"";

            for(const auto &sn : sample_name) {
                os << ",\"" << sn << "\"";
            }

            os << "\n";

            double h = 1.0/(n_samples - 1);

            std::size_t n_vars = sample_name.size();
            for(std::size_t i = 0; i < n_samples; ++i) {
                os << (i*h);

                for(int d = 0; d <  LIBMESH_DIM; ++d) {
                    os << "," << points[i](d);
                }

                const auto &s = sample[i];

                for(std::size_t k = 0; k < n_vars; ++k) {
                    os << "," << s.values[k];
                }

                os << std::endl;
            }

        }

        void export_values() const override
        {
            std::ofstream os(output_path.c_str());

            if(os.good()) {
                describe(os);
                os.close();
            }
        }

        Impl()
        : n_samples(50), output_path("pol.csv") {}

        class Sample {
        public:
            std::vector<double> values;
            explicit Sample(const std::size_t n = 1)
            : values(n, 0)
            {
                assert(n > 0);
            }
        };

    private:
        libMesh::Point p0, p1;
        std::size_t n_samples;

        std::vector<libMesh::Point> points;
        std::vector<std::string> sample_name;
        std::vector<Sample> sample;
        std::string output_path;


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

    template class LinePostProcessor<LibMeshFunctionSpace, UVector>;
}
