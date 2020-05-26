#ifndef UTOPIA_MONITOR_HPP
#define UTOPIA_MONITOR_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    template <class Tensor, int Backend = Traits<Tensor>::Backend>
    class EvalMonitor {
    public:
        using SizeType = typename Traits<Tensor>::SizeType;

        template <class Derived, int Order>
        inline void apply(const SizeType & /*it*/, Tensor & /*t*/
        ) {
            // doing nothing
        }

        template <class Derived, int Order>
        inline void apply(const SizeType & /*it*/,
                          Tensor & /*t*/,
                          const std::string & /*name_of_file*/,
                          const std::string & /*name_of_instance*/
        ) {
            // doing nothing
        }
    };

    /**
     * @defgroup   profiling Profiling/Debugging
     * @ingroup    base_functions
     */

    /**
     * @defgroup   interoperability Interoperability
     * @ingroup    base_functions
     */

    /**
     * @ingroup    profiling
     * @brief      Monitoring of the specific tensor over iterations in solvers. \n
     *             Outputs matlab file, with saved values of tensor over iterates. \n
     *
     * @param[in]  it      The iteration.
     * @param      t       Tensor to be monitored.
     */
    template <class Derived, int Order>
    inline void monitor(const typename Traits<Derived>::SizeType &it, Tensor<Derived, Order> &t) {
        EvalMonitor<Derived>::apply(it, t.derived());
    }

    /**
     * @ingroup    profiling
     * @brief      Monitoring of the specific tensor over iterations in solvers. \n
     *             Outputs matlab file, with saved values of tensor over iterates. \n
     *
     * @param[in]  it                    The iteration number.
     * @param[in]  name_of_file          The name of file.
     * @param[in]  name_of_instance      The name of tensor to be saved.
     * @param      t                     Tensor to be monitored.
     */
    template <class Derived, int Order>
    inline void monitor(const typename Traits<Derived>::SizeType &it,
                        Tensor<Derived, Order> &t,
                        const std::string &name_of_file,
                        const std::string &name_of_instance) {
        EvalMonitor<Derived>::apply(it, t.derived(), name_of_file, name_of_instance);
    }

    /**
     * @ingroup    profiling
     * @brief      Gets number of global nnz.
     *
     * @return     The global nnz.
     */
    template <class Derived>
    inline typename Traits<Derived>::SizeType global_nnz(const Tensor<Derived, 2> &t) {
        return t.global_nnz();
    }

    /**
     * @ingroup    profiling
     * @brief      Gets number of local nnz.
     *
     * @return     The local nnz.
     */
    template <class Derived>
    inline typename Traits<Derived>::SizeType local_nnz(const Tensor<Derived, 2> &t) {
        return t.local_nnz();
    }

}  // namespace utopia

#endif  // UTOPIA_MONITOR_HPP
