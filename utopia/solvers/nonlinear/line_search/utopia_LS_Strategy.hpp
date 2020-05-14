#ifndef UTOPIA_LS_STRATEGY_HPP
#define UTOPIA_LS_STRATEGY_HPP
#include "utopia_FunctionNormalEq.hpp"
#include "utopia_Input.hpp"

namespace  utopia
{
    /**
     * @brief      Base class for different line-search strategies.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */

    template<class Vector>
    class LSStrategy : virtual public Configurable

    {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

    public:
        ~LSStrategy() override {}

        LSStrategy():   verbose_(false),
                        c1_(1e-4),
                        max_it_(50),
                        rho_(0.75),
                        alpha_min_(1e-9)
        {

        };

        /**
         * @brief      Gets the alpha = step-size. Function needs to be provided by actual LS strategies.
         *
         * @param
         * @param[in]
         *
         * @return     The alpha.
         */
        virtual bool get_alpha(FunctionBase<Vector> &, const Vector &, const Vector& , const Vector &, Scalar &) = 0;
        virtual bool get_alpha(LeastSquaresFunctionBase<Vector> &, const Vector &, const Vector& , const Vector &, Scalar &) = 0;


        virtual void verbose(const bool & verbose)
        {
            verbose_ = verbose;
        }


        virtual bool verbose() const
        {
            return verbose_;
        }



        virtual void c1(const Scalar  & c1_in)
        {
            c1_ = c1_in;
        }


        virtual Scalar c1() const
        {
            return c1_;
        }


        virtual void rho(const Scalar  & rho_in)
        {
            rho_ = rho_in;
        }


        virtual Scalar rho() const
        {
            return rho_;
        }

        virtual void alpha_min(const Scalar  & alpha_min)
        {
            alpha_min_ = alpha_min;
        }


        virtual Scalar alpha_min() const
        {
            return alpha_min_;
        }

        virtual void max_it(const SizeType  & max_it)
        {
            max_it_ = max_it;
        }


        virtual SizeType max_it() const
        {
            return max_it_;
        }

        void read(Input &in) override {
            in.get("verbose", verbose_);
            in.get("c1", c1_);
            in.get("max_it", max_it_);
            in.get("rho", rho_);
            in.get("alpha_min", alpha_min_);
        }

        void print_usage(std::ostream &os) const override {
            this->print_param_usage(os, "verbose", "bool", "Verbose.", "false");
            this->print_param_usage(os, "c1", "double", "Constant used for Wolfe conditions.", "1e-4");
            this->print_param_usage(os, "max_it", "int", "Maximum number of iterations.", "50");
            this->print_param_usage(os, "rho", "double", "Contraction factor.", "0.75");
            this->print_param_usage(os, "alpha_min", "double", "Minimum allowed step-size.", "1e-9");
        }

        virtual void init_memory(const Layout &layout) = 0;

    private:
        bool verbose_;      /*!< Verbose inside of LS strategy.  */
        Scalar c1_;         /*!< Constant for Wolfe conditions \f$ c_1 \in (0,1),   c_1 = 10^{-4} \f$.  */
        SizeType max_it_;     /*!< Maximum of the iterations inside of LS strategy.  */
        Scalar rho_;        /*!< Contraction factor.   */
        Scalar alpha_min_;  /*!< Minimum allowed step-size.   */

    };
}

#endif //UTOPIA_LS_STRATEGY_HPP
