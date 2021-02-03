#ifndef UTOPIA_UI_SCALAR_SAMPLER_HPP
#define UTOPIA_UI_SCALAR_SAMPLER_HPP

#include "utopia_fe_base.hpp"

#include "utopia_CSV.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_ui.hpp"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>

#include <libmesh/exact_solution.h>

namespace utopia {

    template <typename Scalar>
    class UIFunction {
    public:
        virtual ~UIFunction() {}
        virtual Scalar eval(const std::vector<Scalar> &x) const = 0;
        virtual bool set_current_block(const int subdomain_id) {
            UTOPIA_UNUSED(subdomain_id);
            return true;
        }
    };

    template <>
    class UIFunction<USerialMatrix> {
    public:
        using Scalar = Traits<USerialMatrix>::Scalar;

        virtual ~UIFunction() {}
        virtual USerialMatrix eval(const std::vector<Scalar> &x) const = 0;
        virtual bool set_current_block(const int subdomain_id) {
            UTOPIA_UNUSED(subdomain_id);
            return true;
        }
    };

    template <typename Scalar>
    class UIConstantFunction final : public UIFunction<Scalar> {
    public:
        virtual ~UIConstantFunction() {}

        UIConstantFunction(const Scalar val) : val_(val) {}

        inline Scalar eval(const std::vector<Scalar> &) const { return val_; }

        const Scalar &value() const { return val_; }

    private:
        Scalar val_;
    };

    template <typename Scalar>
    class UIBoxedFunction final : public UIFunction<Scalar> {
    public:
        UIBoxedFunction(const std::vector<Scalar> &lowbo,
                        const std::vector<Scalar> &upbo,
                        const std::shared_ptr<UIFunction<Scalar>> &fun)
            : lowbo_(lowbo), upbo_(upbo), fun_(fun) {}

        inline Scalar eval(const std::vector<Scalar> &x) const override {
            const auto n = std::min(x.size(), lowbo_.size());
            for (std::size_t i = 0; i < n; ++i) {
                if (lowbo_[i] > x[i] || upbo_[i] < x[i]) {
                    return 0.;
                }
            }

            auto value = fun_->eval(x);
            return value;
        }

    private:
        std::vector<Scalar> lowbo_;
        std::vector<Scalar> upbo_;
        std::shared_ptr<UIFunction<Scalar>> fun_;
    };

    template <typename Scalar>
    class UISubdomainFunction final : public UIFunction<Scalar>, public Configurable {
    public:
        UISubdomainFunction() {}

        inline bool has_default() const { return static_cast<bool>(default_fun_); }

        void set_default(const std::shared_ptr<UIFunction<Scalar>> &fun) { default_fun_ = fun; }

        void read(Input &in) override {
            in.get_all([this](Input &sub_is) {
                std::string type = "constant";
                std::string expr = "0.";
                int block_id = -1;

                Scalar val = 0.0;
                sub_is.get("value", val);
                // sub_is.get("value",  expr);
                sub_is.get("type", type);
                sub_is.get("block", block_id);

                bool fun_is_constant = true;

#ifdef UTOPIA_WITH_TINY_EXPR
                fun_is_constant = type == "constant" || type.empty();

                if (!fun_is_constant) {
                    std::cerr << "Not supported yet" << std::endl;
                    // std::string expr = "x";
                    // is.get("function", expr);
                    // fun = std::make_shared<SymbolicFunction>(expr);
                    // double expr = 1.;
                    // is.get("function", expr);
                    // fun_is_constant = true;
                    // fun = std::make_shared<ConstantCoefficient<double, 0>>(expr);
                }
#else
                assert(fun_is_constant);
#endif  // UTOPIA_WITH_TINY_EXPR

                // const Scalar val = std::atof(expr.c_str());

                if (fun_is_constant) {
                    if (block_id == -1) {
                        default_fun_ = std::make_shared<UIConstantFunction<Scalar>>(val);
                    } else {
                        fun_[block_id] = std::make_shared<UIConstantFunction<Scalar>>(val);
                    }
                }

                std::cout << "value: " << val << " type " << type << " block " << block_id << std::endl;
            });
        }

        inline Scalar eval(const std::vector<Scalar> &x) const override {
            assert(active_fun_);
            auto value = active_fun_->eval(x);
            return value;
        }

        bool set_current_block(const int subdomain_id) override {
            auto it = fun_.find(subdomain_id);
            if (it == fun_.end()) {
                active_fun_ = default_fun_;
                return false;
            } else {
                active_fun_ = it->second;
                return true;
            }
        }

        inline bool has_active_function() const { return static_cast<bool>(active_fun_); }

        inline bool good() const { return static_cast<bool>(default_fun_) || !fun_.empty(); }

    private:
        std::shared_ptr<UIFunction<Scalar>> default_fun_;
        std::shared_ptr<UIFunction<Scalar>> active_fun_;
        std::map<int, std::shared_ptr<UIFunction<Scalar>>> fun_;
    };

    template <typename Scalar>
    class UILambdaFunction final : public UIFunction<Scalar> {
    public:
        using F = std::function<Scalar(const std::vector<Scalar> &)>;

        UILambdaFunction(F fun) : fun_(fun) {}

        inline Scalar eval(const std::vector<Scalar> &x) const override { return fun_(x); }

    private:
        F fun_;
    };

    template <typename F>
    std::shared_ptr<UIFunction<double>> lambda_fun(F fun) {
        return std::make_shared<UILambdaFunction<double>>(fun);
    }

    template <typename F>
    std::shared_ptr<UIFunction<double>> boxed_fun(const std::vector<double> &lowbo,
                                                  const std::vector<double> &upbo,
                                                  F fun) {
        return std::make_shared<UIBoxedFunction<double>>(lowbo, upbo, lambda_fun(fun));
    }

    template <typename Scalar_>
    class ContextFunction<std::vector<Scalar_>, UIFunction<Scalar_>>
        : public Expression<ContextFunction<std::vector<Scalar_>, UIFunction<Scalar_>>> {
    public:
        static const int Order = 0;
        typedef Scalar_ Scalar;

        ContextFunction(const std::shared_ptr<UIFunction<Scalar>> &fun) : fun_(fun) {}

        template <int Backend>
        auto eval(const AssemblyContext<Backend> &ctx) const -> std::vector<Scalar> {
            fun_->set_current_block(ctx.block_id());

            const auto &pts = ctx.fe()[0]->get_xyz();

            const auto n = pts.size();
            std::vector<Scalar> ret(n);

            for (std::size_t i = 0; i < n; ++i) {
                std::vector<Scalar> p = {pts[i](0), pts[i](1), pts[i](2)};
                ret[i] = fun_->eval(p);
            }

            return ret;
        }

    private:
        std::shared_ptr<UIFunction<Scalar>> fun_;
    };

    template <typename T>
    inline ContextFunction<std::vector<T>, UIFunction<T>> ctx_fun(const std::shared_ptr<UIFunction<T>> &fun) {
        return ContextFunction<std::vector<T>, UIFunction<T>>(fun);
    }

    template <typename Scalar>
    class UIScalarSampler final : public Configurable, public UIFunction<Scalar> {
    public:
        UIScalarSampler() {}

        ~UIScalarSampler() {}

        void read(Input &is) override {
            std::string file = "";
            is.get("file", file);

            std::fill(std::begin(min_), std::end(min_), 0.);
            std::fill(std::begin(max_), std::end(max_), 0.);
            std::fill(std::begin(dims_), std::end(dims_), 0);

            is.get("min-x", min_[0]);
            is.get("min-y", min_[1]);
            is.get("min-z", min_[2]);

            is.get("max-x", max_[0]);
            is.get("max-y", max_[1]);
            is.get("max-z", max_[2]);

            is.get("nx", dims_[0]);
            is.get("ny", dims_[1]);
            is.get("nz", dims_[2]);

            n = 0;

            std::size_t n_values = 1;
            for (std::size_t i = 0; i < 3; ++i) {
                range_[i] = max_[i] - min_[i];

                if (dims_[i] != 0) {
                    n++;
                    n_values *= dims_[i];
                }
            }

            values_.reserve(n_values);
            std::ifstream ifs(file.c_str());

            while (ifs.good()) {
                Scalar val;
                ifs >> val;
                values_.push_back(val);
            }

            ifs.close();

            std::size_t n_read_values = values_.size();
            assert(n_read_values == n_values);

            if (n_read_values != n_values) {
                std::cout << "number of values is not consistent"
                             "with the dimensions provided expected "
                          << n_values << " but found " << values_.size() << std::endl;
            }
        }

        Scalar eval(const std::vector<Scalar> &x) const override {
            long ind = index(x);
            if (ind < 0 || ind >= values_.size()) {
                assert(false);
                return 0.;
            }
            // std::cout << ind << " = " << (values_[ind]) << std::endl;
            return values_[ind];
        }

        long index(const std::vector<Scalar> &x) const {
            long ret = 0;

            std::size_t offset = 1;
            for (std::size_t i = 0; i < n; ++i) {
                long ind = std::min(long(dims_[i] - 1), long(floor((x[i] - min_[i]) / range_[i] * dims_[i])));

                if (ind < 0 || ind >= dims_[i]) {
                    // out of range
                    return -1;
                }

                ret += offset * ind;
                offset *= dims_[i];
            }

            assert(std::size_t(ret) < values_.size());
            return ret;
        }

        inline bool empty() const { return values_.empty(); }

        void describe(std::ostream &os) const {
            for (std::size_t i = 0; i < 3; ++i) {
                os << "[" << min_[i] << " " << max_[i] << "] " << range_[i] << " \n";
                os << "dims(" << i << ") " << dims_[i] << "\n";
            }
        }

    private:
        Scalar min_[3];
        Scalar max_[3];
        Scalar range_[3];
        SizeType dims_[3];
        SizeType n;
        std::vector<Scalar> values_;
    };

    template <typename Scalar>
    class Normal final {};

    template <typename Scalar_>
    class ContextFunction<std::vector<USerialVector>, Normal<Scalar_>> final
        : public Expression<ContextFunction<std::vector<USerialVector>, Normal<Scalar_>>> {
    public:
        static const int Order = 1;
        typedef Scalar_ Scalar;

        ContextFunction() {}

        template <int Backend>
        auto eval(const AssemblyContext<Backend> &ctx) const -> std::vector<USerialVector> {
            const auto &n = ctx.fe()[0]->get_normals();
            auto nn = n.size();
            auto dim = ctx.spatial_dimension();

            std::vector<USerialVector> normals(nn);
            for (std::size_t i = 0; i < nn; ++i) {
                normals[i].resize(dim);
                for (int d = 0; d < dim; ++d) {
                    normals[i].set(d, n[i](d));
                }
            }

            return normals;
        }
    };

    inline ContextFunction<std::vector<USerialVector>, Normal<double>> normal() {
        return ContextFunction<std::vector<USerialVector>, Normal<double>>();
    }

    template <>
    class ContextFunction<std::vector<USerialMatrix>, UIFunction<USerialMatrix>>
        : public Expression<ContextFunction<std::vector<USerialMatrix>, UIFunction<USerialMatrix>>> {
    public:
        static const int Order = 2;
        using Scalar = Traits<USerialMatrix>::Scalar;

        ContextFunction(const std::shared_ptr<UIFunction<USerialMatrix>> &fun) : fun_(fun) {}

        template <int Backend>
        auto eval(const AssemblyContext<Backend> &ctx) const -> std::vector<USerialMatrix> {
            fun_->set_current_block(ctx.block_id());

            const auto &pts = ctx.fe()[0]->get_xyz();

            const auto n = pts.size();
            std::vector<USerialMatrix> ret(n);

            for (std::size_t i = 0; i < n; ++i) {
                std::vector<Scalar> p = {pts[i](0), pts[i](1), pts[i](2)};
                ret[i] = fun_->eval(p);
            }

            return ret;
        }

    private:
        std::shared_ptr<UIFunction<USerialMatrix>> fun_;
    };

}  // namespace utopia

#endif  // UTOPIA_UI_SCALAR_SAMPLER_HPP
