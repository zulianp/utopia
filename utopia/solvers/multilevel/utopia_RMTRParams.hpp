#ifndef UTOPIA_RMTR_PARAMS_HPP
#define UTOPIA_RMTR_PARAMS_HPP

#include "utopia_Core.hpp"
#include "utopia_Input.hpp"
#include "utopia_TRBase.hpp"

namespace utopia {

    template <class Vector>
    class RMTRParams : public TrustRegionParams<Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        RMTRParams()
            : _it_global(0),

              _max_coarse_it(2),
              _max_sucessful_smoothing_it(1),
              _max_sucessful_coarse_it(1),
              _max_QP_smoothing_it(5),
              _max_QP_coarse_it(50),
              _eps_delta_termination(0.001),
              _grad_smoothess_termination(1. / 8.),
              _hessian_update_delta(0.15),
              _hessian_update_eta(0.5),
              _use_line_search(false),

              red_(FG_LIGHT_MAGENTA),
              def_(FG_DEFAULT),
              yellow_(FG_LIGHT_YELLOW),
              green_(FG_LIGHT_GREEN) {}

        ~RMTRParams() override = default;

        void read(Input &in) override {
            TrustRegionParams<Vector>::read(in);

            in.get("skip_BC_checks", _skip_BC_checks);
            in.get("max_coarse_it", _max_coarse_it);
            in.get("max_sucessful_smoothing_it", _max_sucessful_smoothing_it);
            in.get("max_sucessful_coarse_it", _max_sucessful_coarse_it);
            in.get("max_QP_smoothing_it", _max_QP_smoothing_it);
            in.get("max_QP_coarse_it", _max_QP_coarse_it);

            in.get("eps_delta_termination", _eps_delta_termination);
            in.get("grad_smoothess_termination", _grad_smoothess_termination);
            in.get("check_gradient_smoothness", _check_gradient_smoothness);

            in.get("use_line_search", _use_line_search);

            in.get("hessian_update_delta", _hessian_update_delta);
            in.get("hessian_update_eta", _hessian_update_eta);

            in.get("deltaH_lagging", _deltaH_lagging);
            in.get("Hpost_lagging", _Hpost_lagging);

            // in.get("deltag_lagging", _deltag_lagging);
        }

        void print_usage(std::ostream &os) const override {
            TrustRegionParams<Vector>::print_usage(os);

            this->print_param_usage(os, "skip_BC_checks", "bool", "Skip treatment of BC conditions.", "false");
            this->print_param_usage(os, "max_coarse_it", "int", "Maximum number of coarse iterations.", "2");
            this->print_param_usage(os, "max_sucessful_smoothing_it", "int", "Maximum number of smoothing steps.", "1");
            this->print_param_usage(os, "max_sucessful_coarse_it", "int", "Maximum number of corse iterations.", "1");

            this->print_param_usage(
                os, "max_QP_smoothing_it", "int", "Maximum number of iterations for fine level QP solver.", "5");
            this->print_param_usage(
                os, "max_QP_coarse_it", "int", "Maximum number of iterations for coarse level QP solver.", "50");

            this->print_param_usage(os,
                                    "eps_delta_termination",
                                    "real",
                                    "Constant used for quiting recursion based on size of tr. radius.",
                                    "0.001");
            this->print_param_usage(os,
                                    "grad_smoothess_termination",
                                    "real",
                                    "Constant used for quiting recursion based on "
                                    "smoothness of the gradient on given level.",
                                    "0.5");
            this->print_param_usage(
                os, "check_gradient_smoothness", "bool", "Flag turning on/off check for gradient smoothiness.", "true");

            this->print_param_usage(
                os, "use_line_search", "bool", "Use line-search for coarse grid correction..", "true");

            this->print_param_usage(os,
                                    "hessian_update_delta",
                                    "real",
                                    "Constant used for deciding whether to assemble fresh hessian or no.",
                                    "0.15");
            this->print_param_usage(os,
                                    "hessian_update_eta",
                                    "real",
                                    "Constant used for deciding whether to assemble fresh hessian or no.",
                                    "0.5");
            this->print_param_usage(
                os, "verbosity_level", "VerbosityLevel", "Specifies level of verbosity.", "VERBOSITY_LEVEL_NORMAL");

            this->print_param_usage(os,
                                    "deltaH_lagging",
                                    "deltaH_lagging",
                                    "Use lagging of delta H term. Applies only for "
                                    "Second order consistent models.",
                                    "False");
            this->print_param_usage(
                os, "Hpost_lagging", "Hpost_lagging", "Use lagging of Hessian for post-smoothing.", "False");

            // this->print_param_usage(os,
            //                         "deltag_lagging",
            //                         "deltag_lagging",
            //                         "Use lagging of gradient for delta g term computation.",
            //                         "False");
        }

        virtual bool deltaH_lagging() const { return _deltaH_lagging; }

        virtual void deltaH_lagging(const bool &deltaH_lagging) { _deltaH_lagging = deltaH_lagging; }

        virtual bool Hpost_lagging() const { return _Hpost_lagging; }

        virtual void Hpost_lagging(const bool &flg) { _Hpost_lagging = flg; }

        // virtual bool deltag_lagging() const { return _deltag_lagging; }

        // virtual void deltag_lagging(const bool &flg) { _deltag_lagging = flg; }

        virtual VerbosityLevel verbosity_level() const { return _verbosity_level; }

        virtual void verbosity_level(const VerbosityLevel &level) { _verbosity_level = level; }

        virtual MultilevelNormSchedule norm_schedule() const { return _norm_schedule; }

        virtual void norm_schedule(const MultilevelNormSchedule &schedule) { _norm_schedule = schedule; }

        virtual void grad_smoothess_termination(const Scalar &grad_smoothess_termination) {
            _grad_smoothess_termination = grad_smoothess_termination;
        }

        virtual Scalar grad_smoothess_termination() const { return _grad_smoothess_termination; }

        virtual void skip_BC_checks(const bool skip_BC_checks) { _skip_BC_checks = skip_BC_checks; }

        virtual bool skip_BC_checks() const { return _skip_BC_checks; }

        virtual void max_coarse_it(const SizeType &max_coarse_it) { _max_coarse_it = max_coarse_it; }

        virtual SizeType max_coarse_it() const { return _max_coarse_it; }

        // this parameter define total number of smoothing its - sucessful and
        // unsuccessful while pre_smooting and post_smoothing count just sucessful
        // iterations
        virtual void max_sucessful_smoothing_it(const SizeType &max_sucessful_smoothing_it) {
            _max_sucessful_smoothing_it = max_sucessful_smoothing_it;
        }

        virtual SizeType max_sucessful_smoothing_it() const { return _max_sucessful_smoothing_it; }

        virtual void max_sucessful_coarse_it(const SizeType &max_sucessful_coarse_it) {
            _max_sucessful_coarse_it = max_sucessful_coarse_it;
        }

        virtual SizeType max_sucessful_coarse_it() const { return _max_sucessful_coarse_it; }

        virtual void max_QP_smoothing_it(const SizeType &num_it) { _max_QP_smoothing_it = num_it; }

        virtual SizeType max_QP_smoothing_it() const { return _max_QP_smoothing_it; }

        virtual void max_QP_coarse_it(const SizeType &num_it) { _max_QP_coarse_it = num_it; }

        virtual SizeType max_QP_coarse_it() const { return _max_QP_coarse_it; }

        virtual void check_grad_smoothness(const bool flg) { _check_gradient_smoothness = flg; }

        virtual bool check_grad_smoothness() const { return _check_gradient_smoothness; }

        virtual void eps_delta_termination(const Scalar &eps) { _eps_delta_termination = eps; }

        virtual Scalar eps_delta_termination() const { return _eps_delta_termination; }

        virtual void hessian_update_delta(const Scalar &delta) { _hessian_update_delta = delta; }

        virtual Scalar hessian_update_delta() const { return _hessian_update_delta; }

        virtual void hessian_update_eta(const Scalar &eps) { _hessian_update_eta = eps; }

        virtual Scalar hessian_update_eta() const { return _hessian_update_eta; }

        virtual bool use_line_search() const { return _use_line_search; }

        virtual void use_line_search(const bool &flg) { _use_line_search = flg; }

    protected:
        SizeType _it_global; /** * global iterate counter  */

        bool _skip_BC_checks{false};

        SizeType _max_coarse_it;              /** * maximum iterations on coarse level   */
        SizeType _max_sucessful_smoothing_it; /** * max smoothing iterations  */
        SizeType _max_sucessful_coarse_it;    /** * max smoothing iterations  */

        SizeType _max_QP_smoothing_it;
        SizeType _max_QP_coarse_it;

        Scalar _eps_delta_termination;      /** * maximum delta allowed on coarse level -
                                               makes sure that coarse level corection stays
                                               inside fine level radius  */
        Scalar _grad_smoothess_termination; /** * determines when gradient is not
                                               smooth enough => does pay off to go to
                                               coarse level at all  */
        bool _check_gradient_smoothness{true};

        Scalar _hessian_update_delta; /** * tolerance used for updating hessians */
        Scalar _hessian_update_eta;   /** * tolerance used for updating hessians */

        bool _use_line_search{false};

        VerbosityLevel _verbosity_level{VERBOSITY_LEVEL_NORMAL};
        MultilevelNormSchedule _norm_schedule{ALWAYS};
        // bool _deltag_lagging{false};

        bool _Hpost_lagging{false};
        bool _deltaH_lagging{false};

        ColorModifier red_;
        ColorModifier def_;
        ColorModifier yellow_;
        ColorModifier green_;
    };

}  // namespace utopia

#endif  // UTOPIA_RMTR_PARAMS_HPP