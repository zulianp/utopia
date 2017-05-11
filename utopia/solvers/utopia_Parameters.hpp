/*
* @Author: alenakopanicakova
* @Date:   2016-05-22
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-10
*/

#ifndef UTOPIA_UTOPIA_PARAMETERS_HPP
#define UTOPIA_UTOPIA_PARAMETERS_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <map>

// ------------------------------------------- WORK IN PROGRESS ---------------------------------------------

namespace utopia 
{
    /**
     * @brief      This class keeps track on all parameters, that we have in linear and nonlinear solvers.
     *             It provides default choice of params and routines to set user-defined preferences. 
     *             
     */
    class Parameters 
    {
        typedef double Scalar;
        typedef long SizeType;

    public:
       Parameters()
       {
          verbose_ = true; 
          time_statistics_ = true; 
          convergence_reason_ = 0; 
          max_it_ = 300; 
          num_it_ = 0; 
          
          tol_ = 1e-7;  // to be depreciated.... 
          

          atol_ = 1e-7; 
          rtol_ = 1e-7; 
          stol_ = 1e-7; 
          

          check_diff_ = false; 
          solver_type_ = "NEWTON"; 

        /*----------  Newton  ----------*/          
          alpha_ = 1; 

        /*----------  TR  ----------*/
          trust_region_alg_ = "STEIHAUG_TOINT"; 
          delta_max_ = 1e8; 
          delta0_ = 1e5; 
          gamma1_ = 0.2; 
          gamma2_ = 2.0; 
          eta1_ = 0.1; 
          eta2_ = 0.85; 
          rho_tol_ = 0.01; 
          SteihaugToint_tol_ = 1e-10; 
          eps_ = 2e-12; 

       /*----------  MG  ----------*/
          blocksize_ = 3; 
          smoother_type_ = 0;
          mg_type_ = 1; 
          pre_smoothing_steps_ = 3; 
          post_smoothing_steps_ = 3; 
          omega_ = 0.66; 
          static_time_step_ = true; 
          cycle_type_      = "multiplicative"; 

        /*----------  LS  ----------*/
          line_search_alg_ = "BACKTRACKING"; 
          ls_rho_ = 0.5; 
          alpha_min_ = 1e-7; 
          c1_ = 1e-4; 
          c2_ = 1e-8; 
          n_line_search_iters_ = 50; 
          line_search_inner_verbose_ = false; 

       /*----------  linear solver  ----------*/
          ksp_atol_ = 1e-8; 
          ksp_rtol_ = 1e-8; 
          ksp_dtol_ = 1e-8; 
          ksp_max_it_ = 100; 
          linear_solver_verbose_ = true;

          linear_solver_type_ = "UTOPIA_CG"; 
          preconditioner_type_ = "lu"; 
          preconditioner_factor_mat_solver_package_ = "mumps"; 

       /*----------  domain decomposition  ----------*/
          overlap_ = 0; 
          local_max_it_ = 20;

      /* -------------  monitoring ---------------- */
          log_iterates_       = false; 
          log_system_         = false; 
          log_norms_          = false; 

          stay_quiet();
        }

    void stay_quiet()
    {
      verbose_ = false;
      time_statistics_ = false;
      linear_solver_verbose_ = false;
    }    

   /*----------------------------  nonlinear solver  ------------------------------------*/
    SizeType num_it() const                  { return num_it_; }
    SizeType convergence_reason() const      { return convergence_reason_; }
    SizeType max_it() const                  { return max_it_; } 
    
    Scalar tol() const                       { return tol_; } 
    Scalar atol() const                       { return atol_; } 
    Scalar stol() const                       { return stol_; } 
    Scalar rtol() const                       { return rtol_; } 
    
    bool  verbose() const                    { return verbose_; } 
    bool  time_statistics() const                    { return time_statistics_; } 
    bool  differentiation_control() const    { return check_diff_; } 
    char const * solver_type()  const        { return solver_type_; }  

    bool  log_iterates() const                { return log_iterates_; } 
    bool  log_system() const                  { return log_system_; } 
    bool  log_norms() const                   { return log_norms_; } 


  /*---------------------------------  Newton  -----------------------------------------------*/
    Scalar alpha() const                       { return alpha_; } 
  /*---------------------------------  TR  -----------------------------------------------*/
    Scalar delta_max() const                  { return delta_max_; } 
    Scalar delta0() const                     { return delta0_; } 
    char const *   trust_region_alg() const   { return trust_region_alg_; }  

    Scalar gamma1() const                     { return gamma1_; } 
    Scalar gamma2() const                     { return gamma2_; }
    Scalar eta1() const                       { return eta1_; } 
    Scalar eta2()  const                      { return eta2_; } 
    Scalar rho_tol() const                    { return rho_tol_; } 
    Scalar SteihaugToint_tol() const          { return SteihaugToint_tol_; } 
    Scalar eps() const                        { return eps_; } 

/*---------------------------------  MG   --------------------------------------------------*/
    SizeType  block_size() const              { return blocksize_; }
    char const *   smoother_type() const      { return smoother_type_; }
    SizeType  mg_type()  const                { return mg_type_; }
    char const *  cycle_type()  const          { return cycle_type_; }
    SizeType  pre_smoothing_steps()const      { return pre_smoothing_steps_; }
    SizeType  post_smoothing_steps()const     { return post_smoothing_steps_; }
    Scalar    omega()const                    { return omega_;    } 
    bool      static_time_step() const        { return static_time_step_; }; 


/*---------------------------------  LS   --------------------------------------------------*/  
    Scalar c1()  const                        { return c1_; } 
    Scalar c2()  const                        { return c2_; } 
    Scalar ls_rho()  const                    { return ls_rho_; } 
    Scalar alpha_min()  const                    { return alpha_min_; } 

    char const * line_search_alg() const      { return line_search_alg_; }  
    SizeType  n_line_search_iters()const      { return n_line_search_iters_; }
    bool line_search_inner_verbose()const     { return line_search_inner_verbose_; }
    
/*---------------------------------  linear solver   -----------------------------------------*/  
    bool      linear_solver_verbose() const  { return linear_solver_verbose_; } 
    Scalar    ksp_atol() const               { return ksp_atol_; } 
    Scalar    ksp_rtol()  const              { return ksp_rtol_; } 
    Scalar    ksp_dtol()  const              { return ksp_dtol_; } 
    SizeType  ksp_max_it()  const            { return ksp_max_it_; } 

    bool      precondition() const          { return precondition_; } 
    bool      linear_solver_time_statistics() const   { return linear_solver_time_statistics_; } 

    char const *      lin_solver_type()  const                            { return linear_solver_type_; }  
    char const *      preconditioner_type()  const                        { return preconditioner_type_; }  
    char const *      preconditioner_factor_mat_solver_package() const    { return preconditioner_factor_mat_solver_package_; }  

  /*-------------------------------  domain decomposition  -------------------------------------*/
    SizeType    overlap()   const            { return overlap_; } 
    SizeType    local_max_it()  const        { return local_max_it_; } 



  /*----------------------------  SETTERS  ------------------------------------*/
    void num_it(const SizeType & num_it)                              { num_it_  = num_it; }
    void convergence_reason(const SizeType & convergence_reason)      { convergence_reason_ = convergence_reason; }
    void max_it(const SizeType & max_it)                              { max_it_ = max_it; } 
    
    void tol(const Scalar & tol )                                       { tol_ = tol; } 
    void rtol(const Scalar & rtol )                                     { rtol_ = rtol; } 
    void stol(const Scalar & stol )                                     { stol_ = stol; } 
    void atol(const Scalar & atol_in )                                  { atol_ = atol_in; } 
    

    void verbose(const bool & verbose)                                { verbose_ = verbose; } 
    void time_statistics(const bool & time_statistics)                { time_statistics_ = time_statistics; } 
    void differentiation_control(const bool & check_diff )            { check_diff_ = check_diff; } 
    void solver_type(char const * solver_type )                       { solver_type_ = solver_type; }  

    void log_iterates(const bool & log_iterates)                      { log_iterates_ = log_iterates; } 
    void log_system(const bool & log_system)                          { log_system_ = log_system; } 
    void log_norms(const bool & log_norms)                            { log_norms_ = log_norms; } 

  /*---------------------------------  Newton  -----------------------------------------------*/
    void alpha(const Scalar & alpha )                                       { alpha_ = alpha; } 

  /*---------------------------------  TR  -----------------------------------------------*/
    void delta_max(const Scalar & delta_max)                          { delta_max_ = delta_max; } 
    void delta0(const Scalar & delta0)                                { delta0_ = delta0; } 
    void trust_region_alg(char const * trust_region_alg)              { trust_region_alg_ = trust_region_alg; }  


    void gamma1(const Scalar & gamma1)                                {  gamma1_ = gamma1; } 
    void gamma2(const Scalar & gamma2)                                {  gamma2_ = gamma2; }
    void eta1(const Scalar & eta1)                                    {  eta1_ = eta1; } 
    void eta2(const Scalar & eta2)                                    {  eta2_ = eta2; } 
    void rho_tol(const Scalar & rho_tol)                              { rho_tol_ = rho_tol; } 
    void SteihaugToint_tol(const Scalar & SteihaugToint_tol)          { SteihaugToint_tol_ = SteihaugToint_tol; } 
    void eps(const Scalar & eps)                                      { eps_ = eps; } 

/*---------------------------------  MG   --------------------------------------------------*/
    void  block_size(const SizeType & block_size)                     {  blocksize_ = block_size; }
    void  smoother_type(char const *  smoother_type)                  {  smoother_type_ = smoother_type; }
    void  mg_type(const SizeType & mg_type)                           {  mg_type_ = mg_type; }
    void  pre_smoothing_steps(const SizeType & pre_smoothing_steps)   {  pre_smoothing_steps_ = pre_smoothing_steps; }
    void  post_smoothing_steps(const SizeType & post_smoothing_steps) {  post_smoothing_steps_ = post_smoothing_steps; }
    void  omega(const SizeType & omega)                               {  omega_ = omega; } 
    void  static_time_step(const bool & static_time_step)             {  static_time_step_ = static_time_step; }; 
    void  cycle_type(char const *  cycle_type)                        {  cycle_type_ = cycle_type; }

/*---------------------------------  LS   --------------------------------------------------*/  
    void  c1(const Scalar & c1)                                       {  c1_ = c1; } 
    void  c2(const Scalar & c2)                                       {  c2_ = c2; } 
    void  ls_rho(const Scalar & ls_rho)                               {  ls_rho_ = ls_rho; } 
    void  alpha_min(const Scalar & alpha_min)                               {  alpha_min_ = alpha_min; } 

    void  line_search_alg( char const * line_search_alg)              {  line_search_alg_ = line_search_alg; }  
    void  n_line_search_iters(const SizeType & n_line_search_iters)   {  n_line_search_iters_ = n_line_search_iters; }
    void line_search_inner_verbose(const bool & verbose)              { line_search_inner_verbose_ = verbose; }
    
/*---------------------------------  linear solver   -----------------------------------------*/  
    void  linear_solver_verbose(const bool & linear_solver_verbose)   {  linear_solver_verbose_ = linear_solver_verbose; } 
    void  ksp_atol(const Scalar & ksp_atol)                           {  ksp_atol_ = ksp_atol; } 
    void  ksp_rtol(const Scalar & ksp_rtol)                           {  ksp_rtol_ = ksp_rtol; } 
    void  ksp_dtol(const Scalar & ksp_dtol)                           {  ksp_dtol_ = ksp_dtol; } 
    void  ksp_max_it(const SizeType & ksp_max_it)                     {  ksp_max_it_ = ksp_max_it; } 

    void lin_solver_type(char const * lin_solver_type)                {  linear_solver_type_ = lin_solver_type; }  
    void preconditioner_type(char const * preconditioner_type)        {  preconditioner_type_ = preconditioner_type; }  
    void preconditioner_factor_mat_solver_package(char const * preconditioner_factor_mat_solver_package)     {  preconditioner_factor_mat_solver_package_ = preconditioner_factor_mat_solver_package; }  

    void  precondition(const bool & precondition)                     {  precondition_ = precondition; }
    void linear_solver_time_statistics(const bool & linear_solver_time_statistics)   { linear_solver_time_statistics_ = linear_solver_time_statistics; }

  /*-------------------------------  domain decomposition  -------------------------------------*/
    void overlap(const SizeType & overlap)                            {  overlap_ = overlap; } 
    void local_max_it(const SizeType & local_max_it)                  {  local_max_it_ = local_max_it; } 


  

    protected: 
          bool verbose_; 
          bool time_statistics_; 
          SizeType convergence_reason_; 
          SizeType max_it_; 
          SizeType num_it_; 
          
          Scalar tol_; 
          Scalar atol_; 
          Scalar stol_; 
          Scalar rtol_;

          Scalar alpha_;

          bool check_diff_; 
          char const  * solver_type_; 

          char const  * trust_region_alg_; 
          Scalar delta_max_; 
          Scalar delta0_; 
          Scalar gamma1_; 
          Scalar gamma2_; 
          Scalar eta1_; 
          Scalar eta2_; 
          Scalar rho_tol_; 
          Scalar SteihaugToint_tol_; 
          Scalar eps_; 

          SizeType  blocksize_; 
          char const  * smoother_type_;
          char const  * cycle_type_;
          SizeType  mg_type_; 
          SizeType  pre_smoothing_steps_; 
          SizeType  post_smoothing_steps_; 
          Scalar    omega_; 
          bool      static_time_step_; 

          char const  * line_search_alg_; 
          Scalar ls_rho_; 
          Scalar c1_; 
          Scalar c2_; 
          Scalar alpha_min_;
          SizeType n_line_search_iters_; 
          bool line_search_inner_verbose_; 

          Scalar ksp_atol_; 
          Scalar ksp_rtol_; 
          Scalar ksp_dtol_; 
          SizeType ksp_max_it_; 
          bool linear_solver_verbose_; 

          char const * linear_solver_type_; 
          char const * preconditioner_type_; 
          char const  * preconditioner_factor_mat_solver_package_; 

          SizeType overlap_; 
          SizeType local_max_it_; 

          bool precondition_ = false; 
          bool linear_solver_time_statistics_ = true; 

          bool log_iterates_; 
          bool log_system_; 
          bool log_norms_;

    };

}


#endif // UTOPIA_UTOPIA_PARAMETERS_HPP
