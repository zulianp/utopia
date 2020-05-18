#ifndef UTOPIA_LIBMESH_LAMBDA_FUNCTION_HPP
#define UTOPIA_LIBMESH_LAMBDA_FUNCTION_HPP

#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"

namespace utopia {

    template <typename Output>
    class LibMeshLambdaFunction : public libMesh::FunctionBase<Output> {
    public:
        typedef libMesh::Real Scalar;
        typedef libMesh::Point Point;

        typedef std::function<Scalar(const Point &)> ScalarFunction;
        typedef std::function<void(const Point &, libMesh::DenseVector<Output> &)> VectorFunction;

        LibMeshLambdaFunction(const VectorFunction fun) : vector_fun_(fun) {}

        LibMeshLambdaFunction(const ScalarFunction fun) : scalar_fun_(fun) {}

        void init() override {}
        void clear() override {}

        libMesh::UniquePtr<libMesh::FunctionBase<Output> > clone() const override {
            //#ifdef LIBMESH_HAVE_CXX14_MAKE_UNIQUE
            //			using libMesh::make_unique;
            //			return make_unique< LibMeshLambdaFunction<Output> >(*this);
            //#else
            return libMesh::UniquePtr<libMesh::FunctionBase<Output> >(new LibMeshLambdaFunction<Output>(*this));
            //#endif
        }

        Output operator()(const Point &p, const Scalar time = 0.) override {
            assert(scalar_fun_);
            return scalar_fun_(p);
        }

        void operator()(const Point &p, const Scalar time, libMesh::DenseVector<Output> &output) override {
            if (!vector_fun_) {
                auto val = (*this)(p, time);

                int size = output.size();
                for (int i = 0; i < size; ++i) {
                    output(i) = val;
                }

            } else {
                vector_fun_(p, output);
            }
        }

        inline libMesh::Real component(unsigned int var_cmp, const Point &p, libMesh::Real time) override {
            if (!vector_fun_) {
                return (*this)(p, time);
            } else {
                libMesh::DenseVector<Output> out;
                (*this)(p, time, out);
                return out(var_cmp);
            }
        }

    private:
        VectorFunction vector_fun_;
        ScalarFunction scalar_fun_;
    };

    template <typename Output, std::size_t N>
    class STL2LibMeshLambdaFunction : public libMesh::FunctionBase<Output> {
    public:
        typedef Output Scalar;
        typedef libMesh::Point Point;

        typedef std::array<Scalar, N> ArrayT;
        typedef std::function<Scalar(const ArrayT &)> ScalarFunction;
        typedef std::function<ArrayT(const ArrayT &)> VectorFunction;

        STL2LibMeshLambdaFunction(const VectorFunction fun) : vector_fun_(fun) {}

        STL2LibMeshLambdaFunction(const ScalarFunction fun) : scalar_fun_(fun) {}

        void init() override {}
        void clear() override {}

        libMesh::UniquePtr<libMesh::FunctionBase<Output> > clone() const override {
            //#ifdef LIBMESH_HAVE_CXX14_MAKE_UNIQUE
            //			using libMesh::make_unique;
            //			return make_unique< STL2LibMeshLambdaFunction<Output> >(*this);
            //#else
            return libMesh::UniquePtr<libMesh::FunctionBase<Output> >(new STL2LibMeshLambdaFunction<Output, N>(*this));
            //#endif
        }

        inline Output operator()(const Point &p, const Scalar time = 0.) override {
            ArrayT in;
            for (std::size_t i = 0; i < N; ++i) {
                in[i] = p(i);
            }

            return scalar_fun_(in);
        }

        inline void operator()(const Point &p, const Scalar time, libMesh::DenseVector<Output> &output) override {
            ArrayT in, out;
            for (std::size_t i = 0; i < N; ++i) {
                in[i] = p(i);
            }

            out = vector_fun_(in);

            for (std::size_t i = 0; i < N; ++i) {
                output(i) = out[i];
            }
        }

        inline libMesh::Real component(unsigned int var_cmp, const Point &p, libMesh::Real time) override {
            assert(var_cmp < N);

            ArrayT in, out;
            for (std::size_t i = 0; i < N; ++i) {
                in[i] = p(i);
            }

            out = vector_fun_(in);
            return out[var_cmp];
        }

    private:
        VectorFunction vector_fun_;
        ScalarFunction scalar_fun_;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_LAMBDA_FUNCTION_HPP
