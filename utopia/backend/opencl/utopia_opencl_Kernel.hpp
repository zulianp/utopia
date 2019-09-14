#ifndef UTOPIA_CL_KERNEL_HPP
#define UTOPIA_CL_KERNEL_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Utils.hpp"
#include "utopia_TreeNavigator.hpp"
#include "utopia_TreeProperties.hpp"

#include "utopia_opencl_CodeTemplate.hpp"
#include "utopia_opencl_Symbol.hpp"
#include "utopia_opencl_Env.hpp"
#include "utopia_opencl_Evaluate.hpp"

#include "utopia_opencl_HasMatMatMul.hpp"
#include "utopia_opencl_Traits.hpp"


namespace utopia {
    namespace opencl {

        class KernelCodeGenerator {
        public:
            template<class Derived>
            void visit(const utopia::Expression<Derived> &expr, Env &env, const Options &options)
            {
                pre_visit();

                std::cout << "FAILED" << std::endl;
                std::cout << "------------------------------------------" << std::endl;
                std::cout << tree_format(expr.derived().get_class()) << std::endl;
                std::cout << "------------------------------------------" << std::endl;
                os_ << "[TODO]";
                assert(false);
                post_visit();
            }

            template<class InnerExpr>
            void visit(const utopia::Transposed<InnerExpr> &expr, Env &env, const Options &options)
            {
                pre_visit();
                visit(expr.expr(), env, make_transposed(options));
                post_visit();
            }

            template<class InnerExpr, class Operation>
            void visit(const utopia::Reduce<InnerExpr, Operation> &expr, Env &env, const Options &options)
            {
                pre_visit();
                visit(expr.expr(), env, remove_transpose(options));
                tpl_.map("REDUCE_OP", Symbol<Operation>::fun_str());
                post_visit();
            }

            template<class Tensor, int Order>
            void visit(const utopia::Wrapper<Tensor, Order> &node, Env &env, const Options &options)
            {
                pre_visit();
                int arg_num = env.new_arg(make_var_wrapper(node));
                std::string var_name = "arg_" + std::to_string(arg_num);
                os_ << var_name;

                if(options.transposed()) {
                    os_ << "[index_transposed]";
                } else {
                    os_ << "[index]";
                }

                if(options.is_const()) {
                    tpl_.map("[list]ARG_IN", "__global const Scalar * " + var_name);
                } else {
                    tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name);
                }

                if(options.is_permutation()) {
                    Size s = size(node);

                    int arg_num_rows = env.new_arg(make_var_copy(Number<int>(s.get(0))));
                    int cols = (s.n_dims() > 1)? s.get(1) : 1;
                    int arg_num_cols = env.new_arg(make_var_copy(Number<int>(cols)));

                    std::string var_name_rows = "arg_" + std::to_string(arg_num_rows);
                    std::string var_name_cols = "arg_" + std::to_string(arg_num_cols);

                    tpl_.map("ROWS", var_name_rows);
                    tpl_.map("COLUMNS", var_name_cols);

                    tpl_.map("[list]ARG_IN", "const int " + var_name_rows);
                    tpl_.map("[list]ARG_IN", "const int " + var_name_cols);
                }

                post_visit();
            }

            template<class Expr>
            void visit(const utopia::Evaluate<Expr> &eval, Env &env, const Options &options)
            {
                pre_visit();
                int arg_num = env.new_arg(make_var_copy(eval));
                std::string var_name = "arg_" + std::to_string(arg_num);
                os_ << var_name;

                if(Expr::Order == 0) {
                    //No
                    tpl_.map("[list]ARG_IN", "const Scalar " + var_name);

                } else {
                    if(Expr::Order == 2 && options.transposed()) {
                        os_ << "[index_transposed]";
                    } else {
                        os_ << "[index]";
                    }

                    if(options.is_const()) {
                        tpl_.map("[list]ARG_IN", "__global const Scalar * " + var_name);
                    } else {
                        tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name);
                    }
                }

                post_visit();
            }

            template<typename T>
            void visit(const utopia::Number<T> &number, Env &env, const Options &options)
            {
                pre_visit();
                int arg_num = env.new_arg(make_var_copy(number));
                std::string var_name = "arg_" + std::to_string(arg_num);
                os_ << var_name;

                tpl_.map("[list]ARG_IN", "const Scalar " + var_name);
                post_visit();
            }

            template<class Left, class Right, class Operation>
            void visit(const utopia::Binary<Left, Right, Operation> &node, Env &env, const Options &options)
            {
                pre_visit();
                os_ << Symbol<Operation>::fun_str();
                os_ << "(";
                visit(node.left(), env, options);
                os_ << ", ";
                visit(node.right(), env, options);
                os_ << ")";
                post_visit();
            }

            template<class Left, class Right>
            void visit(const utopia::Construct<Left, Right> &node, Env &env, const Options &options)
            {
                pre_visit();
                if(has_mat_mat_mul(node.right())) {
                    //TODO
                    visit(node.left(), env, make_permutation(options));
                    os_ << " = (";
                    visit(node.right(), env, make_const(options));
                    os_ << ") ";

                } else {
                    visit(node.left(), env, options);
                    os_ << " = (";
                    visit(node.right(), env, make_const(options));
                    os_ << ") ";
                }
                post_visit();
            }



            template<class Type, int Order>
            void visit(const Factory<Type, Order> &factory, Env &env, const Options &options)
            {
                os_  << "/*Factory method*/ ";
                visit(factory.type(), env, options);
            }

            template<typename T>
            void visit(const Values<T> &value, Env &env, const Options &options)
            {
                int arg_num  = env.new_arg(make_var_copy(utopia::Number<T>(value.value())));
                std::string var_name  = "arg_" + std::to_string(arg_num);
                os_  << var_name;

                tpl_.map("[list]ARG_IN", "const Scalar " + var_name);
            }

            template<int Order>
            void visit(const Factory<Identity, Order> &expr, Env &env, const Options &options)
            {
                Size s = size(expr);
                int arg_num_row  = env.new_arg(make_var_copy(utopia::Number<int>(s.get(0))));
                int arg_num_col  = env.new_arg(make_var_copy(utopia::Number<int>(s.get(1))));

                std::string var_name_row  = "arg_" + std::to_string(arg_num_row);
                std::string var_name_col  = "arg_" + std::to_string(arg_num_col);



                os_  << "identity(" << var_name_row << ", " << var_name_col << ", index)";

                tpl_.map("[list]FUNCTIONS", "Scalar identity(const SizeType rows, const SizeType cols, const SizeType index);\n//");
                tpl_.map("[list]FUNCTIONS", "\nScalar identity(const SizeType rows, const SizeType cols, const SizeType index)\n"
                    "{\n\tconst Scalar ret = (index/cols) == (index % rows) ? 1 : 0;\n"
                    " \treturn ret;\n}");


                tpl_.map("[list]ARG_IN", "const int " + var_name_row);
                tpl_.map("[list]ARG_IN", "const int " + var_name_col);
            }

            template<class Left, class Right>
            void visit(const utopia::Multiply<Left, Right> &node, Env &env, const Options &options)
            {
                pre_visit();

                Size s = size(node.left());

                int arg_num_left  = env.new_arg(make_var_wrapper(node.left()));
                int arg_num_right = env.new_arg(make_var_wrapper(node.right()));
                int arg_num_rows  = env.new_arg(make_var_copy(Number<int>(s.get(0))));
                int arg_num_cols  = env.new_arg(make_var_copy(Number<int>(s.n_dims() == 1? 1 : s.get(1))));

                std::string var_name_left  = "arg_" + std::to_string(arg_num_left);
                std::string var_name_right = "arg_" + std::to_string(arg_num_right);
                std::string var_name_rows  = "arg_" + std::to_string(arg_num_rows);
                std::string var_name_cols  = "arg_" + std::to_string(arg_num_cols);

                if(utopia::order(node) == 2) {
                    Size s = size(node.right());
                    int arg_num_right_cols = env.new_arg(make_var_copy(Number<int>(s.get(1))));
                    std::string var_name_right_cols = "arg_" + std::to_string(arg_num_right_cols);

                    const bool left_transposed  = is_transposed(node.left());
                    const bool right_transposed = is_transposed(node.right());

                    if(left_transposed && right_transposed) {
                        os_ << "matrix_matrix_multiplication_at_entry(";

                        if(options.transposed()) {
                            os_ << "index1, index0, ";
                        } else {
                            os_ << "index0, index1, ";
                        }

                        os_ << var_name_right_cols << ", ";
                        os_ << var_name_cols  << ", ";
                        os_ << var_name_rows  << ", ";
                        os_ << var_name_right << ", ";
                        os_ << var_name_left  << ")";

                    } else {
                        if(left_transposed) {
                            os_ << "matrix_tranposed_matrix_multiplication_at_entry(";
                        } else if(right_transposed) {
                            os_ << "matrix_matrix_transposed_multiplication_at_entry(";
                        } else {
                            os_ << "matrix_matrix_multiplication_at_entry(";
                        }

                        if(options.transposed()) {
                            os_ << "index1, index0, ";
                        } else {
                            os_ << "index0, index1, ";
                        }

                        os_ << var_name_rows << ", ";
                        os_ << var_name_cols  << ", ";
                        os_ << var_name_right_cols  << ", ";
                        os_ << var_name_left << ", ";
                        os_ << var_name_right  << ")";

                        if(options.is_const()) {
                            tpl_.map("[list]ARG_IN", "__global const Scalar * " + var_name_left);
                            tpl_.map("[list]ARG_IN", "__global const Scalar * " + var_name_right);
                        } else {
                            tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name_left);
                            tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name_right);
                        }

                        tpl_.map("[list]ARG_IN", "const int " + var_name_rows);
                        tpl_.map("[list]ARG_IN", "const int " + var_name_cols);
                        tpl_.map("[list]ARG_IN", "const int " + var_name_right_cols);
                    }

                } else {
                    if(options.transposed()) {
                        os_ << "matrix_vector_multiplication_transposed_at_row(index,";
                    } else {
                        os_ << "matrix_vector_multiplication_at_row(index,";
                    }

                    os_ << var_name_rows << ", " << var_name_cols << ", ";
                    os_ << var_name_left << ", " << var_name_right << ")";

                    if(options.is_const()) {
                        tpl_.map("[list]ARG_IN", "__global const Scalar * " + var_name_left);
                        tpl_.map("[list]ARG_IN", "__global const Scalar * " + var_name_right);
                    } else {
                        tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name_left);
                        tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name_right);
                    }

                    // tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name_left);
                    // tpl_.map("[list]ARG_IN", "__global Scalar * " + var_name_right);
                    tpl_.map("[list]ARG_IN", "const int " + var_name_rows);
                    tpl_.map("[list]ARG_IN", "const int " + var_name_cols);
                }

                post_visit();
            }

            template<class InnerExpr, class Operation>
            void visit(const utopia::Unary<InnerExpr, Operation> &node, Env &env, const Options &options)
            {
                pre_visit();
                os_ << Symbol<Operation>::str();
                os_ << "(";
                visit(node.expr(), env, options);
                os_ << ")";
                post_visit();
            }

            template<class Derived>
            bool generate(const Expression<Derived> &expr,
                const std::string &kernel_name,
                const std::string &code_template,
                std::string &code)
            {

                CLStats::instance().code_generation_begin();
                depth_ = 0;
                tpl_.map("KERNEL_NAME", kernel_name);

                visit(expr.derived(), env_, Options());

                // std::cout << tree_format(expr.get_class()) << std::endl;
                // env_.describe(std::cout);


                tpl_.map("OP", os_.str());
                os_.clear();
                const bool ok = tpl_.parse(code_template, code);
                CLStats::instance().code_generation_end();
                return ok;
            }

        private:
            CodeTemplate tpl_;
            std::stringstream os_;
            Env env_;
            int depth_;
            std::string indent_;

            void pre_visit()
            {
                depth_++;
            }

            void post_visit()
            {
                depth_--;
            }
        };

        class KernelNameGenerator {
        public:
            constexpr static const char * base_name()
            {
                return "utopia_auto_generated_kernel_";
            }

            std::string next_name()
            {
                std::string name(base_name());
                return name + std::to_string(n_generated_kernels++);
            }

            static KernelNameGenerator &instance()
            {
                static KernelNameGenerator instance;
                return instance;
            }

        private:
            KernelNameGenerator()
            : n_generated_kernels(0)
            {}

            unsigned long n_generated_kernels;
        };

        class IKernel {
        public:
            virtual ~IKernel() {}
            //remove me
            virtual const std::string &get_code_string() const = 0;
            virtual const std::string &get_name() const = 0;
            virtual bool is_callable() const = 0;
            virtual bool make_callable(cl::Program &program) = 0;
        };




        template<class Expr, int IsScalar = is_scalar_expression_tree<Expr>() >
        class Kernel : public IKernel {
        public:
            static Kernel &instance()
            {
                static Kernel instance;
                return instance;
            }

            bool initialized() const
            {
                return initialized_;
            }

            class RemapVariables {
            public:
                template<class Node>
                void pre_order_visit(const Node &) {}

                template<class Node>
                void in_order_visit(const Node &) {}

                template<class Node>
                void post_order_visit(const Node &) {}

                template<class Tensor, int Order>
                void pre_order_visit(const utopia::Wrapper<Tensor, Order> &node)
                {
                    env_->new_arg(make_var_wrapper(node));
                }

                template<class Left, class Right>
                void in_order_visit(const utopia::Construct<Left, Right> &expr)
                {
                    if(has_mat_mat_mul(expr.right())) {
                        Size s = size(expr.right());
                        for(SizeType i = 0; i < s.n_dims(); ++i) {
                            env_->new_arg(make_var_copy(utopia::Number<int>(s.get(i))));
                        }
                    }
                }

                template<typename T>
                void pre_order_visit(const utopia::Number<T> &node)
                {
                    env_->new_arg(make_var_copy(node));
                }

                template<typename InnerExpr>
                void pre_order_visit(const utopia::Evaluate<InnerExpr> &node)
                {
                    // std::cout << "encountered evaluate for " << node.expr().get_class() << std::endl;
                    env_->new_arg(make_var_copy(node));
                }

                template<class Type, int Order>
                void pre_order_visit(const Factory<Type, Order> &factory)
                {
                    pre_order_visit(factory.type());
                }

                template<typename T>
                void pre_order_visit(const Values<T> &value)
                {
                    env_->new_arg(make_var_copy(utopia::Number<T>(value.value())));
                }

                template<int Order>
                void pre_order_visit(const Factory<Identity, Order> &factory)
                {
                    Size s = size(factory);
                    assert(Order == s.n_dims() && "must have consitend sizes and order");
                    for(SizeType i = 0; i < Order; ++i) {
                        env_->new_arg(make_var_copy(utopia::Number<int>(s.get(i))));
                    }
                }

                template<class Left, class Right>
                void post_order_visit(const Multiply<Left, Right> &expr)
                {
                    Size s = size(expr.left());
                    env_->new_arg(make_var_copy(utopia::Number<int>(s.get(0))));
                    env_->new_arg(make_var_copy(utopia::Number<int>(s.n_dims() == 1? 1 : s.get(1))));

                    if(utopia::order(expr) == 2) {
                        Size s_right = size(expr.right());
                        env_->new_arg(make_var_copy(utopia::Number<int>(s_right.get(1))));
                    }
                }


                RemapVariables(Env &env)
                : env_(&env)
                {
                    env_->clear();
                }

                Env *env_;
            };

            template<class Left, class Right>
            inline std::string get_template_path(const Multiply<Left, Right> &)
            {
                return "/utopia_TemplateBLAS_2.cl";
            }

            template<class InnerExpr, class Operation>
            std::string get_template_path(const Reduce<InnerExpr, Operation> &)
            {
                return "/utopia_TemplateReduce.cl";
            }

            template<class ExprT>
            std::string get_template_path(const ExprT &expr)
            {
                // std::cout << tree_format(expr.get_class()) << std::endl;
                // std::cout << "has_mat_mat_mul: " << has_mat_mat_mul(expr) << std::endl;

                if(has_mat_mat_mul(expr)) {
                    return "/utopia_TemplateBLAS_2.cl";
                } else {
                    return "/utopia_TemplateBLAS_1.cl";
                }
            }

            template<typename T, class InnerExpr>
            void initialize(const Construct< Number<T>,
                Evaluate<InnerExpr> > &)
            {
                static_assert(InnerExpr::Order != 0, "this method should never be called");
                return;
            }

            void initialize(const Expr &expr) {
                if(initialized()) {
                    // std::cout << "initialized" << std::endl;
                    return;
                }

                // if(has_mat_mat_mul(expr)) {
                // 	std::cout << "Needs to handle mat mat mul" << std::endl;
                // }

                KernelCodeGenerator gen;

                std::string code_template;
                bool ok = read(utopia::Utopia::instance().get("opencl_templates_path") + get_template_path(expr), code_template); assert(ok);

                name_ = KernelNameGenerator::instance().next_name();

                ok = gen.generate(expr, name_, code_template, code_); assert(ok);
                // std::cout << "initializing kernel for:\n" << treeFormat(expr.get_class()) << std::endl;
                initialized_ = true;
            }

            inline std::pair<cl_ulong, cl_ulong> optimal_group_setup(const cl_ulong n)
            {
                // cl_ulong wgSize = 512;
                cl_ulong preferred_work_group_size;
                kernel_.getWorkGroupInfo(CLContext::instance().current_device(), CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, &preferred_work_group_size);
                // clGetKernelWorkGroupInfo(kernel_, CLContext::instance().current_device(), CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(cl_ulong), &preferred_work_group_size, NULL);
                cl_ulong group_size = std::min(n, preferred_work_group_size);
                cl_ulong n_groups   = std::max(n/group_size, cl_ulong(1));
                return std::make_pair(n_groups, group_size);
            }

            template<typename T, class InnerExpr, class Operation>
            bool call(const Evaluate<
                Construct<
                Number<T>,
                Evaluate< Reduce<InnerExpr, Operation> >
                >
                > &eval)
            {


                const auto &constr = eval.expr();

                if(!call(constr.right())) {
                    return false;
                }

                constr.left() = constr.right().get_value();
                return true;
            }


            template<class InnerExpr, class Operation>
            bool call(const Evaluate< Reduce<InnerExpr, Operation> > &eval)
            {

                CLStats::instance().kernel_execution_begin();

                const auto &expr = eval.expr();

                RemapVariables remap(env_);
                auto nav = make_nav(remap);
                nav.set_prune_evaluate_branch(true);
                // nav.setVerbose(true);
                nav.visit(expr);

                typedef UTOPIA_SCALAR(InnerExpr) Scalar;

                cl_int err = CL_SUCCESS;

                auto &queue = CLContext::instance().current_queue();

                Size s = size(expr.expr());

                int n_entries_to_reduce = s.get(0);
                for(SizeType i = 1; i < s.n_dims(); ++i) {
                    n_entries_to_reduce *= s.get(i);
                }

                auto setup = optimal_group_setup(n_entries_to_reduce);
                cl_ulong n_groups   = setup.first;
                cl_ulong group_size = setup.second;


                assert(n_entries_to_reduce > 0);
                assert(n_groups > 0);
                CLVector<Scalar> result(n_groups);

                err = kernel_.setArg(0, n_entries_to_reduce); 							assert(check_cl_error(err));
                err = kernel_.setArg(1, result.buffer); 								assert(check_cl_error(err));
                err = kernel_.setArg(2, cl::__local(sizeof(Scalar) * group_size)); 		assert(check_cl_error(err));

                // std::cout << "additional args: " << env_.args().size() << std::endl;

                for(const auto &arg_ptr : env_.args()) {
                    arg_ptr->set_as_arg_at(arg_ptr->get_arg_num() + 3, kernel_);		assert(check_cl_error(err));
                }

                err = queue.enqueueNDRangeKernel(kernel_, cl::NullRange,
                    cl::NDRange(n_entries_to_reduce),
                    cl::NDRange(group_size)); 											assert(check_cl_error(err));
                err = queue.finish(); 													assert(check_cl_error(err));

                result.synch_read_buffer(queue);

                Scalar host_reduced = result.entries[0];

                for(SizeType i = 1; i < result.entries.size(); ++i) {
                    host_reduced = Operation::template apply<Scalar>(host_reduced, result.entries[i]);
                }

                eval.set_value(host_reduced);

                // std::cout << "[REDUCTION] = " << host_reduced << std::endl;

                CLStats::instance().kernel_execution_end();
                return true;
            }



            // template<class ExprT>
            // bool call(const Evaluate<ExprT> &eval)
            // {
            // 	auto w = wrap<ExprT::Order>(eval.backend_tensor_ptr());
            // 	return call( construct(w, eval.expr()) );

            // }

            template<class Left, class Right>
            bool call(const Evaluate< Construct<Left, Right> > &eval)
            {
                return call(eval.expr());
            }

            inline cl::NDRange opt_global_range(const Size &s) const
            {

                auto &dev = CLContext::instance().current_device();

                if(CL_DEVICE_TYPE_GPU == dev.getInfo<CL_DEVICE_TYPE>()) {
                    int d_2 = (s.n_dims() == 2)? s.get(1) : 1;
                    return cl::NDRange(s.get(0), d_2, 1);
                } else {
                    int max_compute_units = dev.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
                    int g_work_size_0 = s.get(0);
                    g_work_size_0 = std::min(max_compute_units, g_work_size_0);
                    int g_work_size_1 = 1;

                    if(s.n_dims() >= 2) {
                        g_work_size_1 = s.get(1);
                        g_work_size_1 = std::min(std::max(1, max_compute_units - g_work_size_0), g_work_size_1);
                    }

                    return cl::NDRange(g_work_size_0, g_work_size_1, 1);
                }
            }


            template<class Derived>
            bool call(const Expression<Derived> &all_expr)
            {
                CLStats::instance().kernel_execution_begin();

                cl_int err = CL_SUCCESS;

                const Derived &expr = all_expr.derived();

                RemapVariables remap(env_);
                auto nav = make_nav(remap);
                nav.set_prune_evaluate_branch(true);
                // nav.setVerbose(true);
                nav.visit(expr);

                Size s = size(expr);
                // std::cout << "expr size" << std::endl;
                // disp(s);

                // std::cout << "additional args: " << env_.args().size() << std::endl;

                int n_entries = s.get(0);
                for(SizeType i = 1; i < s.n_dims(); ++i) {
                    n_entries *= s.get(i);
                }

                cl::NDRange global_range = opt_global_range(s);

                err = kernel_.setArg(0, n_entries); assert(check_cl_error(err));

                assert(n_entries > 0);
                // std::cout << "n_args: " << (env_.args().size() + 1) << std::endl;

                for(const auto &arg_ptr : env_.args()) {
                    // std::cout << "setting arg: " << arg_ptr->get_arg_num() << std::endl;
                    // err = kernel_.setArg(arg_ptr->get_arg_num() + 1, arg_ptr->buffer());
                    arg_ptr->set_as_arg_at(arg_ptr->get_arg_num() + 1, kernel_);
                    assert(check_cl_error(err));
                }

                // std::cout << "evaluating: " << get_name() << std::endl;
                err = CLContext::instance().current_queue().enqueueNDRangeKernel(kernel_, cl::NullRange, global_range, cl::NullRange); assert(check_cl_error(err));
                err = CLContext::instance().current_queue().finish(); assert(check_cl_error(err));

                CLStats::instance().kernel_execution_end();
                return true;
            }

            bool execute(const Expr &expr)
            {
                if(initialized()) {
                    // std::cout << "reusing kernel remap variables" << std::endl;
                    RemapVariables remap(env_);
                    make_nav(remap).visit(expr);
                } else {
                    initialize(expr);
                }


                // cl::CommandQueue queue(CLContext::instance().current(), CLContext::instance().current_device());

                // //write arrays A and B to the device
                // queue.enqueueWriteBuffer(buffer_A,CL_TRUE,0,sizeof(int)*10,A);
                // queue.enqueueWriteBuffer(buffer_B,CL_TRUE,0,sizeof(int)*10,B);

                // //run the kernel
                // cl::KernelFunctor callback(cl::Kernel(program, get_name().c_str()), queue, cl::NullRange, cl::NDRange(10), cl::NullRange);
                // callback(buffer_A,buffer_B,buffer_C);


                //alternative way to run the kernel
                /*
                       cl::Kernel kernel_add=cl::Kernel(program,"simple_add");
                       kernel_add.setArg(0,buffer_A);
                       kernel_add.setArg(1,buffer_B);
                       kernel_add.setArg(2,buffer_C);
                       queue.enqueueNDRangeKernel(kernel_add,cl::NullRange,cl::NDRange(10),cl::NullRange);
                       queue.finish();
                 */
                    }

                    const std::string &get_code_string() const override
                    {
                        return code_;
                    }

                    const std::string &get_name() const override
                    {
                        return name_;
                    }

                    bool is_callable() const override
                    {
                        return is_callable_;
                    }

                    bool make_callable(cl::Program &program) override
                    {
                        assert(!is_callable() && "do not call this if it is already callable");

                        kernel_ = cl::Kernel(program, get_name().c_str());
                        is_callable_ = true;
                        return is_callable_;
                    }

                private:
                    Kernel() : initialized_(false), is_callable_(false) {}
                    bool initialized_;
                    bool is_callable_;
                    std::string name_;
                    std::string code_;
                    Env env_;

                    cl::Kernel kernel_;
                };



        template<class Expr>
                class Kernel<Expr, 1> : public IKernel {
                public:
                    typedef typename Expr::Scalar Scalar;

                    void initialize(const Expr &) {}

                    class EvalAction {
                    public:
                template<class Node>
                        void pre_order_visit(const Node &) {

                        }

                template<class Node>
                        void post_order_visit(const Node &) {

                        }

                template<class InnerExpr>
                        void post_order_visit(const Evaluate<InnerExpr> &expr)
                        {
                            value_ = expr.get_value();
                        }

                template<class InnerExpr, class Operation>
                        void post_order_visit(const Unary<InnerExpr, Operation> &expr)
                        {
                            value_ = Operation::template Fun<Scalar>()(value_);
                        }

                template<class Tensor>
                        void post_order_visit(const Wrapper<Tensor, 0> &expr)
                        {
                            value_ = expr.get_value();
                        }

                template<class Node>
                        void in_order_visit(const Node &) {
                    //store the left value for binary type expressions
                            left_value_ = value_;
                        }

                template<class Left, class Right>
                        void post_order_visit(const Multiply<Left, Right> &)
                        {
                            value_ *= left_value_;
                        }

                template<class Left, class Right, class Operation>
                        void post_order_visit(const Binary<Left, Right, Operation> &)
                        {
                            value_ = Operation::template Fun<Scalar>()(left_value_, value_);
                        }

                        Scalar get_value() const
                        {
                            return value_;
                        }

                        EvalAction() : value_(0) {}

                    private:
                        Scalar left_value_;
                        Scalar value_;
                    };

            template<class EvaluatedExpr>
                    bool call(const Evaluate< Construct<Number<Scalar>,
                        Evaluate<EvaluatedExpr>
                        > > &expr)
                    {
                        expr.expr().left() = expr.expr().right().get_value();
                        return true;
                    }



            template<class EvaluatedExpr>
                    bool call(const Evaluate<EvaluatedExpr> &expr)
                    {
                // std::cout << tree_format(expr.get_class()) << std::endl;

                        EvalAction eval;
                        auto nav = make_nav(eval);
                        nav.set_prune_evaluate_branch(true);
                        nav.visit(expr.expr());
                        expr.set_value(eval.get_value());
                        return true;
                    }


                    bool call(const Expr &expr)
                    {
                // std::cout << "HERE: is scalar" << std::endl;
                        return false;
                    }

                    static Kernel &instance()
                    {
                        static Kernel instance;
                        return instance;
                    }

                    const std::string &get_code_string() const
                    {
                        return null_;
                    }

                    const std::string &get_name() const
                    {
                        return null_;
                    }


                    bool is_callable() const
                    {
                        return false;
                    }

                    bool make_callable(cl::Program &)
                    {
                        return false;
                    }

                private:
                    std::string  null_;
                    Kernel() {}
                };

        template<class Expr>
                Kernel<Expr> &get_kernel(const Expr &expr)
                {
                    Kernel<Expr>::instance().initialize(expr);
                    return Kernel<Expr>::instance();
                }

        template<class Expr>
                bool execute_kernel(const Expr &expr)
                {
                    return Kernel<Expr>::instance().execute(expr);
                }
            }
        }


#endif //UTOPIA_CL_KERNEL_HPP
