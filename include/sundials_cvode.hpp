#ifndef SUNDIALS_CVODE_STRATEGY_HPP
#define SUNDIALS_CVODE_STRATEGY_HPP

#include <iostream>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>
#include <functional>
#include <memory>

#include <cassert>

#include <cvode/cvode.h>                  // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>       // access to serial N_Vector  
#include <sunmatrix/sunmatrix_dense.h>    // access to dense SUNMatrix 
#include <sunlinsol/sunlinsol_dense.h>    // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_band.h>     // access to band SUNMatrix
#include <sunlinsol/sunlinsol_band.h>     // access to band SUNLinearSolver
#include <sunlinsol/sunlinsol_spgmr.h>    // access to SPGMR SUNLinearSolver
#include <sunlinsol/sunlinsol_spfgmr.h>   // access to SPGMR SUNLinearSolver
#include <sunlinsol/sunlinsol_spbcgs.h>   // access to SPBCGS SUNLinearSolver
#include <sunlinsol/sunlinsol_sptfqmr.h>  // access to SPTFQMR SUNLinearSolver 
#include <sundials/sundials_types.h>      // defs. of realtype, sunindextype

#include <fmt/format.h>

#include "sundials_types.hpp"

namespace SUNDIALS {

  /* 
     ---------------------------------------------------------------------------
     CVODE inetrafce
     Implementation: strategy pattern
     Strategies:
     - Direct methods:
       * Dense matrix
       * Band matrix
     - Iterative methods:
       * SPGMR   = Scaled Preconditioned Generalized Minimal Residual
       * SPFGMR  = Scaled Preconditioned Flexible Generalized Minimal
       * SPBCGS  = Scaled Preconditioned Bi-Conjugate Gradient Stable
       * SPTFQMR = Scaled Preconditioned Transpose-Free Quasi-Minimal Residual 
     ---------------------------------------------------------------------------
  */
  
  // ------------ CVODE strategy interface
  class ICVODEStrategy {

  protected:
    
    void* memory;
    sunindextype n_states;
    N_Vector state;
    SUNLinearSolver linear_solver;
    
  private:
    
    Types::CommonUserFunctions common_user_functions;
    Types::LinearMultisptepMethod linear_multi_step_meth;
    N_Vector constraints;
    Types::IntegrationTime time;
    Types::IntegrationTolerance tolerance;
    void* user_data;
    bool is_initialised;
    class IsSet {
    public:
      bool initial_conditions;
      bool constraints;
      bool integration_time;
      bool common_user_functions;
      IsSet()
        : initial_conditions{ false }
        , constraints{ false }
        , integration_time{ false }
        , common_user_functions{ false }
      { }
    } is_set;
    
    virtual int SetLinearSolver() = 0;
    virtual int BindUserFunctions() = 0;
    
  public:

    ICVODEStrategy()
      : memory{ nullptr }
      , n_states{ -1 }
      , state{ nullptr }
      , linear_solver{ nullptr }
      , common_user_functions()
      , linear_multi_step_meth{ Types::LinearMultisptepMethod::invalid }
      , constraints{ nullptr }
      , time()
      , tolerance()
      , user_data{ nullptr }
      , is_initialised{ false }
      , is_set()
    {}

    virtual ~ICVODEStrategy() {
      if( state ) N_VDestroy_Serial( state );
      if( constraints ) N_VDestroy_Serial( constraints );
      if( memory ) CVodeFree( &memory );
      if( linear_solver ) SUNLinSolFree( linear_solver );
    }

    void SetLinearMultiStepMethod( Types::LinearMultisptepMethod const linear_multi_step_meth_ ) {
      linear_multi_step_meth = linear_multi_step_meth_;
    }

    void SetUserFucntions( CVRhsFn right_hand_side
                           , CVProjFn projection
                           , CVErrHandlerFn error_handler ) {
      assert( right_hand_side );
      common_user_functions.right_hand_side = right_hand_side;
      common_user_functions.projection      = projection;
      common_user_functions.error_handler   = error_handler;
      is_set.common_user_functions = true;
    }

    void SetUserData( void* const user_data_ ) {
      user_data = user_data_;
    }

    void SetInitialConditions( std::vector<realtype> const& initial_conditions) {
      assert( !initial_conditions.empty() );
      n_states = initial_conditions.size();
      state = N_VNew_Serial( n_states );
      assert( state );
      realtype* const state_ptr = N_VGetArrayPointer( state );
      std::copy( initial_conditions.begin()
                 , initial_conditions.end()
                 , state_ptr );
      is_set.initial_conditions = true;
    }

    void SetConstraints( std::vector<Types::Constraint> const& constraints_vec ) {
      assert( is_set.initial_conditions );
      assert( n_states = constraints_vec.size() );
      constraints = N_VNew_Serial( n_states );
      assert( state );
      std::vector<double> con_trans;
      std::transform( constraints_vec.begin()
                      , constraints_vec.end()
                      , std::back_inserter( con_trans )
                      , []( Types::Constraint const con )
                      { return static_cast<double>( con ); } );
      realtype* const constraints_ptr = N_VGetArrayPointer( constraints );
      std::copy( con_trans.begin()
                 , con_trans.end()
                 , constraints_ptr );
      is_set.constraints = true;
    }
  
    void SetTolerances( realtype const relative_tolerance
                        , realtype const absolute_tolerance ) {
      tolerance.relative = relative_tolerance;
      tolerance.absolute = absolute_tolerance;
    }

    void SetIntegrationTime( realtype const time_start
                             , realtype const time_stop ) {
      assert( time_start >= 0.0 );
      assert( time_stop > time_start );
      time.start = time_start;
      time.stop = time_stop;
      is_set.integration_time = true;
    }


    /*
      TODO
      Setters aee needed for: 
        CVodeSetMaxOrd
        CVodeSetMaxNumSteps
        CVodeSetMaxHnilWarns
        CVodeSetStabLimDet
        CVodeSetInitStep
        CVodeSetMinStep
        CVodeSetMaxStep
        CVodeSetStopTime
        CVodeSetMaxErrTestFails
        CVodeSetMaxNonlinIters
        CVodeSetMaxConvFails
        CVodeSetNonlinConvCoef
    */
    
    int Initialise() {
      int ret_code;

       // instantiate the solver
      memory = CVodeCreate( static_cast<int>( linear_multi_step_meth ) );
      assert( memory );
      
      // set user data
      if( user_data ) {
        ret_code = CVodeSetUserData( memory
                                     , user_data );
        if( ret_code != CV_SUCCESS ) return ret_code;
      }

      // set the error handler
      if( common_user_functions.error_handler ) {
        ret_code = CVodeSetErrHandlerFn( memory
                                         , common_user_functions.error_handler
                                         , nullptr );
        if( ret_code != CV_SUCCESS ) return ret_code;
      }

      // initialize the integrator memory and specify the user's
      // right hand side function in u'= f(t,u), the inital time,
      // and the initial dependent variable vector u
      ret_code = CVodeInit( memory
                            , common_user_functions.right_hand_side
                            , time.start
                            , state );
      if( ret_code != 0 ) return ret_code;

      // set the linear solver
      ret_code = SetLinearSolver();
      if( ret_code != CV_SUCCESS ) return ret_code;

      // set UDFs here
      ret_code = BindUserFunctions();
      if( ret_code != CV_SUCCESS ) return ret_code;
      
      // set the projection function
      if( common_user_functions.projection ) {
        ret_code = CVodeSetProjFn( memory
                                   , common_user_functions.projection );
        if( ret_code != CV_SUCCESS ) return ret_code;
        // hardcode: disable error projection
        ret_code = CVodeSetProjErrEst( memory, 0 );
        if( ret_code != CV_SUCCESS ) return ret_code;
        // hardcode: projection frequency
        ret_code = CVodeSetProjFrequency(memory, 1);
        if( ret_code != CV_SUCCESS ) return ret_code;
        ret_code =  CVodeSetEpsProj(memory, 1.0e-5);
        if( ret_code != CV_SUCCESS ) return ret_code;
      }
    
      // set the constraints if any
      if( is_set.constraints ) {
        ret_code = CVodeSetConstraints( memory, constraints );
        if( ret_code != CV_SUCCESS ) return ret_code;
      }
      
      
      // set the tolerances
      ret_code = CVodeSStolerances( memory
                                    , tolerance.relative
                                    , tolerance.absolute );
      if( ret_code != CV_SUCCESS ) return ret_code;

      // set the stop time
      ret_code = CVodeSetStopTime( memory, time.stop );
      if( ret_code != CV_SUCCESS ) return ret_code;
      
      //ret_code = CVodeSetMaxNumSteps( memory, 300000 );
      // if( ret_code != CV_SUCCESS ) return ret_code;
      //ret_code = CVodeSetMaxStep( memory, 1.0e-5);
      // if( ret_code != CV_SUCCESS ) return ret_code;

      is_initialised = true;

      return 0;
      
    }

    bool Initialised() const {
      return is_initialised;
    }

    int ReInitialise( double const time_init,
                      std::vector<double> const& state_init_vec ) {
      assert( is_initialised );
      assert( n_states == state_init_vec.size() );
      realtype * const state_init_ptr = N_VGetArrayPointer( state );
      std::copy( state_init_vec.begin()
                 , state_init_vec.end()
                 , state_init_ptr );
      return CVodeReInit( memory
                          , time_init
                          , state );
    }
  
    // TODO: something better than copying twice needs to be done
    int Integrate( realtype const time_target
                   , std::vector<realtype>& state_vec ) {
      assert( is_initialised );
      assert( n_states == state_vec.size() );
      realtype * const state_ptr = N_VGetArrayPointer( state );
      std::copy( state_vec.begin()
                 , state_vec.end()
                 , state_ptr );
      int const ret_code = CVode( memory
                                  , time_target
                                  , state
                                  , &time.reached
                                  , CV_NORMAL ); 
      std::copy( state_ptr
                 , state_ptr + n_states
                 , state_vec.begin() );
      return ret_code;
    }

    int MainSolverStatistics( Types::MainSolverStatistics& stats ) {
      int ret_code;
      assert( memory );
      ret_code = CVodeGetNumSteps( memory, &stats.n_internal_steps );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumRhsEvals( memory, &stats.n_rhs_evals );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumLinSolvSetups( memory, &stats.n_linear_solver_setups );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumErrTestFails( memory, &stats.n_error_test_fails );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumNonlinSolvIters( memory, &stats.n_nonlinear_solver_iters );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumNonlinSolvConvFails( memory, &stats.n_nonlinear_solver_conv_fails );
      if( ret_code != CV_SUCCESS ) return ret_code;      
      ret_code = CVodeGetLastOrder( memory, &stats.last_order_used );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetCurrentOrder( memory, &stats.current_order );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetLastStep( memory, &stats.last_internal_step_size );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetCurrentStep( memory, &stats.next_internal_step_size );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetActualInitStep( memory, &stats.first_internal_step_size );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumStabLimOrderReds( memory, &stats.n_stability_order_reductions );
      if( ret_code != CV_SUCCESS ) return ret_code;
      return 0;
    }

    int LinearSolverStatistics( Types::LinearSolverStatistics& stats ) {
      int ret_code;
      assert( memory );
      ret_code = CVodeGetNumJacEvals( memory, &stats.n_jac_evals );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumLinRhsEvals( memory, &stats.n_rhs_jac_evals );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumLinIters( memory, &stats.n_linear_iterations );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumLinConvFails( memory, &stats.n_linear_conv_fails);
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumPrecEvals( memory, &stats.n_prec_evals );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumPrecSolves( memory, &stats.n_prec_solves );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumJTSetupEvals( memory, &stats.n_jac_vec_setup_evals );
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code = CVodeGetNumJtimesEvals( memory, &stats.n_jac_vec_prod_evals );
      if( ret_code != CV_SUCCESS ) return ret_code;
      return 0;
    }

    std::string PrintSolverStatistics() {
      std::string str;
      {
        Types::MainSolverStatistics stats;
        MainSolverStatistics( stats );
        str += fmt::format( "\nMain solver statistics:\n"
                            "Cumulative number of internal steps         {}\n"
                            "No. calls to r.h.s. function                {}\n"
                            "No. calls to linear solver setup function   {}\n"
                            "No. local error test failures               {}\n"
                            "No. nonlinear solver iterations             {}\n"
                            "No. nonlinear convergence failures          {}\n"
                            "Order used during the last step             {}\n"
                            "Order to be attempted on the next step      {}\n"
                            "Step size used for the last step            {}\n"
                            "Step size to be attempted on the next step  {}\n"
                            "Actual initial step size used               {}\n"
                            "No. order reds due to stab lim detection    {}\n"
                            , stats.n_internal_steps
                            , stats.n_rhs_evals
                            , stats.n_linear_solver_setups
                            , stats.n_error_test_fails
                            , stats.n_nonlinear_solver_iters
                            , stats.n_nonlinear_solver_conv_fails
                            , stats.last_order_used
                            , stats.current_order
                            , stats.last_internal_step_size
                            , stats.next_internal_step_size
                            , stats.first_internal_step_size
                            , stats.n_stability_order_reductions );
      }
      
      {
        Types::LinearSolverStatistics stats;
        LinearSolverStatistics( stats );
        str += fmt::format( "\nLinear solver statistics:\n"
                            "No. Jac evaluations {}\n"
                            "No. RHS. calls for Jac[-vec] evals      {}\n"
                            "No. linear iterations                   {}\n"
                            "No. linear convergence failures         {}\n"
                            "No. preconditioner evaluations          {}\n"
                            "No. preconditioner solves               {}\n"
                            "No. Jacobian-vector setup evaluations   {}\n"
                            "No. Jacobian-vector product evaluations {}\n"
                            , stats.n_jac_evals
                            , stats.n_rhs_jac_evals
                            , stats.n_linear_iterations
                            , stats.n_linear_conv_fails
                            , stats.n_prec_evals
                            , stats.n_prec_solves
                            , stats.n_jac_vec_setup_evals
                            , stats.n_jac_vec_prod_evals );
      }
      
      return str;
    }

  };

  // ------------ Direct solvers interface
  class Direct: public virtual ICVODEStrategy {
  public:
    Direct()
      : ICVODEStrategy()
      , matrix{ nullptr }
      , user_functions()
    {}
    
    ~Direct() {
      if( matrix ) SUNMatDestroy(matrix);
    }
        
    void SetUserFucntions( CVLsJacFn jacobian ) {
      // jacobian is set by the lib if not provided
      user_functions.jacobian = jacobian;
    }
    
  protected:
     SUNMatrix matrix;

  private:
    Types::DirectUserFunctions user_functions;
    
    int BindUserFunctions( ) override {
      return CVodeSetJacFn( memory
                            , user_functions.jacobian );
    }
    
  };

  // ------------ Dense matrix direct strategy
  class Dense : public virtual ICVODEStrategy
              , public Direct {
  public:
    
    Dense()
      : ICVODEStrategy()
      , Direct()
    {}
    
    ~Dense() {}
  private:
    
    int SetLinearSolver() override {
      // create dense SUNMatrix for use in linear solves
      matrix = SUNDenseMatrix( n_states, n_states );
      // create the linear solver
      linear_solver = SUNLinSol_Dense( state, matrix );
      return CVodeSetLinearSolver( memory
                                   , linear_solver
                                   , matrix );
    }
    
  };

  // ------------ Band matrix direct strategy
  class Band : public virtual ICVODEStrategy
             , public Direct {
  public:
    Band()
      : ICVODEStrategy()
      , Direct()
      , bandwidth()
    {}
    
    ~Band() {}

    void SetMatrixBandwidth( Types::Bandwidth const bandwidth_ ) {
      bandwidth = bandwidth_;
    }
    
    void SetMatrixBandwidth( sunindextype const upper
			     , sunindextype const lower ) {
      bandwidth.upper = upper;
      bandwidth.lower = lower;
    }

  private:
    
    Types::Bandwidth bandwidth;

    int SetLinearSolver() override {
      // create band SUNMatrix for use in linear solves
      assert( bandwidth.upper >= 0 );
      assert( bandwidth.lower >= 0 );
      matrix = SUNBandMatrix( n_states
                              , bandwidth.upper
                              , bandwidth.lower );
      // create the linear solver
      linear_solver = SUNLinSol_Band( state, matrix );
      return CVodeSetLinearSolver( memory
                                   , linear_solver
                                   , matrix );
    }
    
  };
 
  // ------------ Iterative solvers interface
  class Iterative: public virtual ICVODEStrategy {
  
  public:
    Iterative()
      : ICVODEStrategy()
      , user_functions()
      , options()
    { }
    
    ~Iterative() {}
    
    void SetUserFucntions( CVLsJacTimesSetupFn jacobian_setup
                           , CVLsJacTimesVecFn jacobian_times_vector
                           , CVLsPrecSetupFn preconditioner_setup
                           , CVLsPrecSolveFn preconditioner_solve ) {
      user_functions.jacobian_setup        = jacobian_setup;
      user_functions.jacobian_times_vector = jacobian_times_vector;
      assert( preconditioner_setup );
      assert( preconditioner_solve );
      user_functions.preconditioner_setup = preconditioner_setup;
      user_functions.preconditioner_solve = preconditioner_solve;
      //is_set.user_functions = true;
    }
    
    void SetLinearSolverOptions( Types::Preconditioner const preconditioner
                                 , Types::GramSchmidt const gram_schmidt
                                 , int const n_krylov_basis_vectors
                                 , int const max_restarts ) {
      options = { preconditioner
                  , gram_schmidt
                  , n_krylov_basis_vectors
                  , max_restarts };
    }
    
  private:

    Types::IterativeUserFunctions user_functions;
            
    int BindUserFunctions() override {
      // set the Jacobian functions
      int ret_code;

      ret_code = CVodeSetJacTimes( ICVODEStrategy::memory
                                   , user_functions.jacobian_setup
                                   , user_functions.jacobian_times_vector );
      if( ret_code != 0 ) return ret_code;
      
      // Set the preconditioner functions
      ret_code = CVodeSetPreconditioner( memory
                                         , user_functions.preconditioner_setup
                                         , user_functions.preconditioner_solve );
      if( ret_code != 0 ) return ret_code;
      
      return ret_code;
    }

  protected:
    Types::IterativeLinarSolverOptions options;
  };

  // ------------ SPGMR iterative strategy
  class SPGMR : public virtual ICVODEStrategy
              , public Iterative {
  public:
    
    SPGMR( )
      : ICVODEStrategy( )
      , Iterative( )
    { }
    
    ~SPGMR( ) { }
    
  private:
    
    int SetLinearSolver( ) override {
      // create the linear solver
      linear_solver = SUNLinSol_SPGMR( state
                                       , static_cast<int>( options.preconditioner )
                                       , options.n_krylov_basis_vectors );
      if( static_cast<void*>( linear_solver) == nullptr ) return -1;

      // set the Gram-Schmidt orthogonilisation method
      int ret_code = SUNLinSol_SPGMRSetGSType( linear_solver
                                               , static_cast<int>( options.gram_schmidt ) );
      if( ret_code != 0 ) return ret_code;

      // set maximum number of allowable restarts
      ret_code = SUNLinSol_SPGMRSetMaxRestarts( linear_solver
                                                , options.max_restarts );
      if( ret_code != 0 ) return ret_code;
      
      // attach the linear solver to the memory
      return CVodeSetLinearSolver( memory
                                   , linear_solver
                                   , nullptr );
    }
  };

  // ------------ SPFGMR iterative strategy
  class SPFGMR : public virtual ICVODEStrategy
               , public Iterative {
  public:
    
    SPFGMR( )
      : ICVODEStrategy( )
      , Iterative( )
    { }
    
    ~SPFGMR( ) { }
    
  private:
    
    int SetLinearSolver( ) override {
      // create the linear solver
      linear_solver = SUNLinSol_SPFGMR( state
                                       , static_cast<int>( options.preconditioner )
                                       , options.n_krylov_basis_vectors );
      if( static_cast<void*>( linear_solver) == nullptr ) return -1;

      // set the Gram-Schmidt orthogonilisation method
      int ret_code = SUNLinSol_SPFGMRSetGSType( linear_solver
                                               , static_cast<int>( options.gram_schmidt ) );
      if( ret_code != 0 ) return ret_code;

      // set maximum number of allowable restarts
      ret_code = SUNLinSol_SPFGMRSetMaxRestarts( linear_solver
                                                , options.max_restarts );
      if( ret_code != 0 ) return ret_code;
      
      // attach the linear solver to the memory
      return CVodeSetLinearSolver( memory
                                   , linear_solver
                                   , nullptr );
    }
  };

  // ------------ SPBCGS iterative strategy
  class SPBCGS : public virtual ICVODEStrategy
               , public Iterative {
  public:
    
    SPBCGS( )
      : ICVODEStrategy( )
      , Iterative( )
    { }
    
    ~SPBCGS( ) { }
    
  private:
    
    int SetLinearSolver( ) override {
      // create the linear solver
      linear_solver = SUNLinSol_SPBCGS( state
                                       , static_cast<int>( options.preconditioner )
                                       , options.n_krylov_basis_vectors );
      if( static_cast<void*>( linear_solver) == nullptr ) return -1;

      // set the Gram-Schmidt orthogonilisation method
      int ret_code = SUNLinSol_SPBCGSSetPrecType( linear_solver
                                                  , static_cast<int>( options.gram_schmidt ) );
      if( ret_code != 0 ) return ret_code;

      // set maximum number of allowable restarts
      ret_code = SUNLinSol_SPBCGSSetMaxl( linear_solver
                                          , options.max_restarts );
      if( ret_code != 0 ) return ret_code;
      
      // attach the linear solver to the memory
      return CVodeSetLinearSolver( memory
                                   , linear_solver
                                   , nullptr );  
    }
  };
  
  // ------------ SPTFQMR iterative strategy
  class SPTFQMR : public virtual ICVODEStrategy
                , public Iterative {
  public:
    
    SPTFQMR( )
      : ICVODEStrategy( )
      , Iterative( )
    { }
    
    ~SPTFQMR( ) { }
    
  private:
    
    int SetLinearSolver( ) override {

      // create the linear solver
      linear_solver = SUNLinSol_SPTFQMR( state
                                       , static_cast<int>( options.preconditioner )
                                       , options.n_krylov_basis_vectors );
      if( static_cast<void*>( linear_solver) == nullptr ) return -1;

      // no Gram-Schmidt orthogonilisation method
      // set maximum number of allowable restarts
      int ret_code = SUNLinSol_SPTFQMRSetMaxl( linear_solver
                                               , options.max_restarts );
      if( ret_code != 0 ) return ret_code;
      
      // attach the linear solver to the memory
      return CVodeSetLinearSolver( memory
                                   , linear_solver
                                   , nullptr );
    }
  };
  

  // CVODE client - context
  class CVODE {
    
  private:
    
    ICVODEStrategy* strategy;

  public:

    explicit CVODE()
      : strategy( nullptr )
    {}
    
    explicit CVODE( ICVODEStrategy* strategy_ )
      : strategy( strategy_ )
    {}
    
    void SetStrategy( ICVODEStrategy* strategy_ ) {
      strategy = strategy_;
    }

    void SetLinearMultiStepMethod( Types::LinearMultisptepMethod const linear_multi_step_meth ) {
      strategy->SetLinearMultiStepMethod( linear_multi_step_meth );
    }
    
    void SetUserFucntions( CVRhsFn right_hand_side
                           , CVProjFn projection
                           , CVErrHandlerFn error_handler ) {
      strategy->SetUserFucntions( right_hand_side
                                  , projection
                                  , error_handler );
    }

    void SetUserData( void* const user_data ) {
      strategy->SetUserData( user_data );
    }

    void SetInitialConditions( std::vector<realtype> const& initial_conditions) const {
      strategy->SetInitialConditions( initial_conditions );
    }

    void SetConstraints( std::vector<Types::Constraint> const& constraints_vec ) const {
      strategy->SetConstraints( constraints_vec );
    }
  
    void SetTolerances( realtype const relative_tolerance
                        , realtype const absolute_tolerance ) const {
      strategy->SetTolerances( relative_tolerance, absolute_tolerance );
    }

    void SetIntegrationTime( realtype const time_start
                             , realtype const time_stop ) const {
      strategy->SetIntegrationTime( time_start, time_stop );
    }
    
    int Initialise() const {
      return strategy->Initialise();
    }
    
    int Integrate( realtype const time_target
                   , std::vector<realtype>& state_vec ) const {
      return strategy->Integrate( time_target, state_vec );
    }


    int MainSolverStatistics( Types::MainSolverStatistics& stats ) const {
      return strategy->MainSolverStatistics( stats );
    }

    int LinearSolverStatistics( Types::LinearSolverStatistics& stats ) const {
      return strategy->LinearSolverStatistics( stats );
    }
    
    std::string PrintSolverStatistics() const {
      return strategy->PrintSolverStatistics( );
    }
  };

} // namespace SUNDIALS


#endif // SUNDIALS_CVODE_STRATEGY_HPP
