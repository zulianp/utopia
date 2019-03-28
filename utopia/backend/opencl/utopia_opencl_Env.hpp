#ifndef UTOPIA_ENV_HPP
#define UTOPIA_ENV_HPP

#include "utopia_opencl_Traits.hpp"

#include "cl.hpp"
#include <vector>


namespace utopia {
    namespace opencl {

        class Options {
        public:
            Options()
            : _transposed(false), is_permutation_(false), is_const_(false)
            {}

            inline bool transposed() const
            {
                return _transposed;
            }

            inline bool is_permutation() const
            {
                return is_permutation_;
            }

            inline bool is_const() const
            {
                return is_const_;
            }


            inline friend Options make_transposed(const Options &opts)
            {
                Options new_opts = opts;
                new_opts.toggleTransposed();
                return new_opts;
            }

            inline friend Options remove_transpose(const Options &opts)
            {
                Options new_opts = opts;
                new_opts.setTransposed(false);
                return new_opts;
            }

            inline friend Options make_permutation(const Options &opts)
            {
                Options new_opts = opts;
                new_opts.is_permutation_ = true;
                return new_opts;
            }

            inline friend Options remove_permutation(const Options &opts)
            {
                Options new_opts = opts;
                new_opts.is_permutation_ = false;
                return new_opts;
            }

            inline friend Options make_const(const Options &opts)
            {
                Options new_opts = opts;
                new_opts.is_const_ = true;
                return new_opts;
            }

        private:
            bool _transposed;
            bool is_permutation_;
            bool is_const_;

            inline void toggleTransposed()
            {
                _transposed = !_transposed;
            }

            inline void setTransposed(const bool transposed)
            {
                _transposed = transposed;
            }
        };


        class Var {
        public:
            inline void set_arg_num(const int arg_num)
            {
                arg_num_ = arg_num;
            }

            inline int get_arg_num() const
            {
                return arg_num_;
            }

            inline long get_id() const
            {
                return id_;
            }

            inline void set_id(const long id)
            {
                id_ = id;
            }

            inline bool has_null_id() const
            {
                return id_ == -1;
            }

            virtual ~Var() {}

            Var()
            :  arg_num_(-1), id_(-1), aliased_arg_num_(-1)
            {}

            inline bool is_alias() const
            {
                return get_aliased_arg_num() != -1;
            }

            inline int get_aliased_arg_num() const
            {
                return aliased_arg_num_;
            }

            inline void set_aliased_arg_num(const int aliased_arg_num)
            {
                aliased_arg_num_ = aliased_arg_num;
            }

            virtual bool set_as_arg_at(const int arg_num, cl::Kernel &Kernel) = 0;


            // virtual const cl::Memory &buffer() const = 0;

            virtual void describe(std::ostream &os) const = 0;


        private:
            int arg_num_;
            long id_;
            int aliased_arg_num_;
        };

        template<typename T, int Order>
        cl::Buffer &get_buffer(Wrapper<T, Order> &v)
        {
            return raw_type(v).buffer;
        }

        template<typename T, int Order>
        const cl::Buffer &get_buffer(const Wrapper<T, Order> &v)
        {
            return raw_type(v).buffer;
        }

        // template<typename T>
        // const cl::Buffer &get_buffer(const Wrapper<CLVector<T>, 1> &v)
        // {
        // 	return v.implementation().buffer;
        // }

        // template<typename T>
        // cl::Buffer &get_buffer(Wrapper<CLMatrix<T>, 2> &v)
        // {
        // 	return v.implementation().buffer;
        // }

        // template<typename T>
        // const cl::Buffer &get_buffer(const Wrapper<CLMatrix<T>, 2> &v)
        // {
        // 	return v.implementation().buffer;
        // }

        template<class Expr>
        class VarRef : public Var {
        public:
            VarRef(const Expr &expr)
            : expr_(expr)
            {}

            const cl::Buffer &buffer() const
            {
                return get_buffer(expr_);
            }

            virtual bool set_as_arg_at(const int arg_num, cl::Kernel &Kernel)
            {
                const cl_int err = Kernel.setArg(arg_num, buffer());  assert(check_cl_error(err));
                return CL_SUCCESS == err;
            }

            virtual void describe(std::ostream &os) const
            {
                os << expr_.getClass() << "\n";
            }

        private:
            UTOPIA_STORE_CONST(Expr) expr_;
        };

        template<class Expr, int Order = Expr::Order>
        class VarCopy : public Var {
        public:
            VarCopy(const Expr &expr)
            : expr_(expr)
            {}

            const cl::Buffer &buffer() const
            {
                return get_buffer(expr_);
            }

            virtual bool set_as_arg_at(const int arg_num, cl::Kernel &Kernel)
            {
                const cl_int err = Kernel.setArg(arg_num, buffer()); assert(check_cl_error(err));
                return CL_SUCCESS == err;
            }

            virtual void describe(std::ostream &os) const
            {
                os << expr_.getClass() << "\n";
            }

        private:
            Expr expr_;
        };

        template<class Expr>
        class VarCopy<Evaluate<Expr, 0>, 0> : public Var {
        public:
                typedef UTOPIA_SCALAR(Expr) Scalar;

                VarCopy(const Evaluate<Expr, 0> &expr)
                : expr_(expr)
                {}

                virtual bool set_as_arg_at(const int arg_num, cl::Kernel &Kernel)
                {
                    const cl_int err = Kernel.setArg(arg_num, expr_.get_value()); assert(check_cl_error(err));
                    return CL_SUCCESS == err;
                }

                virtual void describe(std::ostream &os) const
                {
                    os << expr_.getClass() << "\n";
                }

            private:
                Evaluate<Expr, 0> expr_;
        };


        template<typename T>
        class VarCopy< utopia::Number<T>, 0> : public Var {
        public:
                VarCopy(const utopia::Number<T> &expr)
                : value_(expr)
                {}


                virtual bool set_as_arg_at(const int arg_num, cl::Kernel &Kernel)
                {
                    const cl_int err = Kernel.setArg(arg_num, value_); assert(check_cl_error(err));
                    return CL_SUCCESS == err;
                }

                virtual void describe(std::ostream &os) const
                {
                    os << value_ << "\n";
                }

            private:
                T value_;

        };

        template<class Expr>
        std::shared_ptr< VarCopy<Expr> > make_var_copy(const Expr &expr)
        {
            return std::make_shared< VarCopy<Expr> >(expr);
        }

        template<class Expr>
        std::shared_ptr< VarRef<Expr> > make_var_wrapper(const Expr &expr)
        {
            return std::make_shared< VarRef<Expr> >(expr);
        }

        class Env {
        public:
            int new_arg(const std::shared_ptr<Var> &arg)
            {
                arg->set_arg_num(args_.size());
                args_.push_back(arg);
                return args_.size();
            }

            void detect_aliases(Var &arg)
            {
                using std::min;

                if(!arg.has_null_id()) {
                    bool found_once = false;
                    int arg_num = -1;

                    for(auto  it = args_.begin(); it != args_.end(); ++it) {
                        auto arg_ptr = *it;

                        if(arg.get_id() == arg_ptr->get_id()) {
                            if(found_once) {
                                arg_ptr->set_aliased_arg_num(arg_num);
                            } else {
                                found_once = true;
                                arg_num = arg_ptr->get_arg_num();
                            }
                        } else {

                        }
                    }
                }
            }

            void clear()
            {
                args_.clear();
            }

            std::vector< std::shared_ptr<Var> > &args()
            {
                return args_;
            }

            void describe(std::ostream &os) const
            {
                os << "------------------------------\n";
                for(const auto arg_ptr : args_) {
                    arg_ptr->describe(os);
                }
                os << "------------------------------\n";
            }

        private:
            std::vector< std::shared_ptr<Var> > args_;
        };
    }
}

#endif //UTOPIA_ENV_HPP
