#include <iostream>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>
#include <functional>
#include <memory>
#include <stdexcept>

#include <cassert>

#include "sundials_types.hpp"
#include "sundials_cvode.hpp"
#include "sundials_client.hpp"

namespace SUNDIALS {
  
  ICVODEStrategy::ICVODEStrategy()
    : memory{ nullptr }
    , client{ nullptr }
    , linear_solver{ nullptr }
    , linear_multi_step_meth{ Types::LinearMultisptepMethod::invalid }
    , time_step_ctrl()
    , solver_ctrl()
    , is_initialised{ false }
    , must_be_reinitialised{ false }
    {}
  
  ICVODEStrategy::~ICVODEStrategy() {
    if( memory ) CVodeFree( &memory );
    if( linear_solver ) SUNLinSolFree( linear_solver );
  }

  void ICVODEStrategy::SetLinearMultiStepMethod( Types::LinearMultisptepMethod
                                                 const linear_multi_step_meth_ ) {
    linear_multi_step_meth = linear_multi_step_meth_;
  }

  void ICVODEStrategy::SetClientData( void* const client_ ) {
    if( !client_ ) {
      throw std::invalid_argument( "ICVODEStrategy::SetClientData: "
                                   "The client data is null.");
    }
    client = static_cast<Client* const>( client_ );
  }

  boost::optional<int> ICVODEStrategy::SetInitStepSize( realtype const init_step ) {
    if( init_step < 0.0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetInitStepSize: "
                                   "The initial step size must be positive.");
    }
    time_step_ctrl.init_step = init_step;
    return ( is_initialised
             ? CVodeSetInitStep( memory, time_step_ctrl.init_step )
             : boost::optional<int>() );
  }

  boost::optional<int> ICVODEStrategy::SetMinStepSize( realtype const min_step ) {
    if( min_step < 0.0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMinStepSize: "
                                   "The minimum step size must be positive.");
    }
    time_step_ctrl.min_step = min_step;
    return ( is_initialised
             ? CVodeSetMinStep( memory, time_step_ctrl.min_step )
             : boost::optional<int>() );
  }
  
  boost::optional<int> ICVODEStrategy::SetMaxStepSize( realtype const max_step ) {
    if( max_step < 0.0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMaxStepSize: "
                                   "The maximum step size must be positive.");
    }
    time_step_ctrl.max_step = max_step;
    return ( is_initialised
             ? CVodeSetMaxStep( memory, time_step_ctrl.max_step )
             : boost::optional<int>() );
  }

  boost::optional<int> ICVODEStrategy::SetMaxNumSteps( long int const n_steps ) {
    if( n_steps < 0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMaxNumSteps: "
                                   "The maximum number of steps must be positive.");
    }
    time_step_ctrl.n_steps = n_steps;
    return ( is_initialised
             ? CVodeSetMaxNumSteps( memory, time_step_ctrl.n_steps )
             : boost::optional<int>() );
  }

  void ICVODEStrategy::SetMaxOrder( int const max_order ) {
    if( max_order < 0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMaxOrder: "
                                   "The maximum order of the linear multistep "
                                   "method must be positive.");
    }
    solver_ctrl.max_order = max_order;
    // once set, the maximum order cannot be reset
  }
  
  boost::optional<int> ICVODEStrategy::SetMaxWarnMessages( int const max_warn_msgs ) {
    if( max_warn_msgs < 0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMaxWarnMessages: "
                                   "The maximum number of warning messages "
                                   "must be positive.");
    }
    solver_ctrl.max_warn_msgs = max_warn_msgs;
    return ( is_initialised
             ? CVodeSetMaxHnilWarns( memory, solver_ctrl.max_warn_msgs )
             : boost::optional<int>() );
  }
  
  boost::optional<int> ICVODEStrategy::SetStabilityLimitDetection( bool const stab_lim_det_active ) {
    solver_ctrl.stab_lim_det_active = stab_lim_det_active;
    if( is_initialised ) {
      if( Types::LinearMultisptepMethod::BDF == linear_multi_step_meth ) {
        if( solver_ctrl.stab_lim_det_active && solver_ctrl.max_order < 3 ) {
          // TODO: log info message stating that this call is recommended when
          //       solver_ctrl.max_order >= 3 where BDF is not A-stable
        }
        return CVodeSetStabLimDet( memory
                                   , static_cast<booleantype>( solver_ctrl.stab_lim_det_active ) );
      } else {
        // TODO: log warning message stating that this call has no effect because
        //       Types::LinearMultisptepMethod::Adams == linear_multi_step_meth
        return boost::optional<int>();
      }
    } else {
      return boost::optional<int>();
    }
  }
  
  boost::optional<int> ICVODEStrategy::SetMaxErrorTestFailures( int const max_err_test_fails ) {
    if( max_err_test_fails <= 0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMaxErrorTestFailures: "
                                   "The maximum number of error test failures "
                                   "must be strictly positive.");
    }
    solver_ctrl.max_err_test_fails = max_err_test_fails;
    return ( is_initialised
             ? CVodeSetMaxErrTestFails( memory, solver_ctrl.max_err_test_fails )
             : boost::optional<int>() );
  }
  
  boost::optional<int> ICVODEStrategy::SetMaxNonlinearIterations( int const max_nonlin_iters ) {
    if( max_nonlin_iters <= 0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMaxNonlinearIterations: "
                                   "The maximum number of nonlinear iterations "
                                   "must be strictly positive.");
    }
    solver_ctrl.max_nonlin_iters = max_nonlin_iters;
    return ( is_initialised
             ? CVodeSetMaxNonlinIters( memory, solver_ctrl.max_nonlin_iters )
             : boost::optional<int>() );
  }
  
  boost::optional<int> ICVODEStrategy::SetMaxConvergenceFailures( int const max_conv_fails ) {
    if( max_conv_fails <= 0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetMaxConvergenceFailures: "
                                   "The maximum number of nonlinear iterations "
                                   "must be strictly positive.");
    }
    solver_ctrl.max_conv_fails = max_conv_fails;
    return ( is_initialised
             ? CVodeSetMaxConvFails( memory, solver_ctrl.max_conv_fails )
             : boost::optional<int>() );
  }
  
  boost::optional<int> ICVODEStrategy::SetNonlinConvergenceCoefficient( realtype const nonlin_conv_coef ) {
    if( nonlin_conv_coef <= 0.0 ) {
      throw std::invalid_argument( "ICVODEStrategy::SetNonlinConvergenceCoefficient: "
                                   "The safety factor used in the nonlinear convergence "
                                   "test must  be strictly positive.");
    }
    solver_ctrl.nonlin_conv_coef = nonlin_conv_coef;
    return ( is_initialised
             ? CVodeSetNonlinConvCoef( memory, solver_ctrl.nonlin_conv_coef )
             : boost::optional<int>() );
  }
  
  int ICVODEStrategy::Initialise() {
    int ret_code;

    // instantiate the solver
    memory = CVodeCreate( static_cast<int>( linear_multi_step_meth ) );
    assert( memory );
      
    // set client data
    assert( client );
    ret_code = CVodeSetUserData( memory
                                 , client );
    if( ret_code != CV_SUCCESS ) return ret_code;
    
    // set the error handler
    if( client->opt_udf_set.error_handler ) {
      ret_code = CVodeSetErrHandlerFn( memory
                                       , Client::ErrorHandlerCallback
                                       , nullptr );
      if( ret_code != CV_SUCCESS ) return ret_code;
    }

    // initialize the integrator memory and specify the user's
    // right hand side function in u'= f(t,u), the inital time,
    // and the initial dependent variable vector u
    ret_code = CVodeInit( memory
                          , Client::RightHandSideCallback
                          , client->time.start
                          , client->state );
    if( ret_code != 0 ) return ret_code;

    // set the linear solver
    ret_code = SetLinearSolver();
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set UDFs here
    ret_code = BindUserFunctions();
    if( ret_code != CV_SUCCESS ) return ret_code;
      
    // set the projection function
    if( client->opt_udf_set.projection ) {
      ret_code = CVodeSetProjFn( memory
                                 , Client::ProjectionCallback );
      if( ret_code != CV_SUCCESS ) return ret_code;
      // hardcoded: disable error projection
      ret_code = CVodeSetProjErrEst( memory, 0 );
      if( ret_code != CV_SUCCESS ) return ret_code;
      // hardcoded: projection frequency
      ret_code = CVodeSetProjFrequency(memory, 1);
      if( ret_code != CV_SUCCESS ) return ret_code;
      ret_code =  CVodeSetEpsProj(memory, 1.0e-5);
      if( ret_code != CV_SUCCESS ) return ret_code;
    }
    
    // set the constraints if any
    if( client->has_constraints ) {
      ret_code = CVodeSetConstraints( memory, client->constraint );
      if( ret_code != CV_SUCCESS ) return ret_code;
    }
    
    // set the tolerances
    if(  1 == client->tolerance.absolute.size() ) {
      ret_code = CVodeSStolerances( memory
                                    , client->tolerance.relative
                                    , client->tolerance.absolute.front() );
    } else {
      N_Vector absolute_tolerance = N_VMake_Serial( client->n_states
                                                    , client->tolerance.absolute.data());
      //N_VNew_Serial( client->n_states );
      ret_code = CVodeSVtolerances( memory
				    , client->tolerance.relative
				    , absolute_tolerance );
      N_VDestroy( absolute_tolerance );
    }
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the stop time
    ret_code = CVodeSetStopTime( memory, client->time.stop );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the initial step size
    ret_code = CVodeSetInitStep( memory, time_step_ctrl.init_step );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the minimum step size
    ret_code = CVodeSetMinStep( memory, time_step_ctrl.min_step );
    if( ret_code != CV_SUCCESS ) return ret_code;
    
    // set the maximum step size
    ret_code = CVodeSetMaxStep( memory, time_step_ctrl.max_step );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the maximum number of internal steps (subdivisions) within a time step 
    ret_code = CVodeSetMaxNumSteps( memory, time_step_ctrl.n_steps );
    if( ret_code != CV_SUCCESS ) return ret_code;
    
    // set max order of the linear multistep method
    if( -1 == solver_ctrl.max_order ) {
      solver_ctrl.max_order =
        ( Types::LinearMultisptepMethod::BDF == linear_multi_step_meth )
        ? 5
        : 12;
    }
    ret_code = CVodeSetMaxOrd( memory, solver_ctrl.max_order );
    if( ret_code != CV_SUCCESS ) return ret_code;
    
    // set the maximum number of warning messages
    ret_code = CVodeSetMaxHnilWarns( memory, solver_ctrl.max_warn_msgs );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the activity of the stability
    ret_code = CVodeSetStabLimDet( memory
                                   , static_cast<booleantype>( solver_ctrl.stab_lim_det_active ) );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the maximum number of error test failures
    ret_code = CVodeSetMaxErrTestFails( memory, solver_ctrl.max_err_test_fails );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the maximum number of nonlinear iteartions
    ret_code = CVodeSetMaxNonlinIters( memory, solver_ctrl.max_nonlin_iters );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the maximum number of nonlinear convergence failures
    ret_code = CVodeSetMaxConvFails( memory, solver_ctrl.max_conv_fails );
    if( ret_code != CV_SUCCESS ) return ret_code;

    // set the safety factor used in the nonlinear convergence test
    ret_code = CVodeSetNonlinConvCoef( memory, solver_ctrl.nonlin_conv_coef );
    if( ret_code != CV_SUCCESS ) return ret_code;   
    
    // last statement in this method
    is_initialised = true;

    return 0;
      
  }

  int ICVODEStrategy::ReInitialise( realtype const time,
                                    std::vector<realtype> const& state ) {
    assert( is_initialised );
    assert( client->n_states == state.size() );
    realtype * const state_ptr = N_VGetArrayPointer( client->state );
    std::copy( state.begin()
               , state.end()
               , state_ptr );
    return CVodeReInit( memory
                        , time
                        , client->state );
  }

  int ICVODEStrategy::ReInitialise( realtype const time,
                                    realtype const * const state ) {
    assert( is_initialised );
    realtype * const state_ptr = N_VGetArrayPointer( client->state );
    std::copy( state
               , state + client->n_states
               , state_ptr );
    return CVodeReInit( memory
                        , time
                        , client->state );
  }
  
  /* 
     This method returns a success, warning or failure code to the client code.
     The client code must:
      - decide what to do in case of failure
      - assess the impact of warnings when they occurr
  */
  int ICVODEStrategy::Integrate( realtype const time_target ) {
    assert( is_initialised );
    int ret_code = CVode( memory
                          , time_target
                          , client->state
                          , &client->time.reached
                          , CV_NORMAL );
    /*
      As a last ditch effort, can include an internal recovery mechanism here
      before passing control to the client code:

      if( CV_SUCCESS != ret_code ) {
        possible types of errors:
         - recoverable
         - unrecoverable
      }
      If internal recovery fails, the client code can still 
      take any necessary measures, re-initialise and attempt to integrate again
     */
    return ret_code;
  }

  realtype* ICVODEStrategy::Solution() const {
    return N_VGetArrayPointer( client->state );
  }

  realtype ICVODEStrategy::Solution( sunindextype const i ) const {
    return N_VGetArrayPointer( client->state )[i];
  }

  int ICVODEStrategy::MainSolverStatistics( Types::MainSolverStatistics& stats ) const {
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

  int ICVODEStrategy::LinearSolverStatistics( Types::LinearSolverStatistics& stats ) const {
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

  std::string ICVODEStrategy::PrintSolverStatistics() const {
    std::string str;
    {
      Types::MainSolverStatistics stats;
      MainSolverStatistics( stats );
      str += fmt::format( "\nMain solver statistics:\n"
                          "Cumulative number of internal steps        {}\n"
                          "No. calls to r.h.s. function               {}\n"
                          "No. calls to linear solver setup function  {}\n"
                          "No. local error test failures              {}\n"
                          "No. nonlinear solver iterations            {}\n"
                          "No. nonlinear convergence failures         {}\n"
                          "Order used during the last step            {}\n"
                          "Order to be attempted on the next step     {}\n"
                          "Step size used for the last step           {}\n"
                          "Step size to be attempted on the next step {}\n"
                          "Actual initial step size used              {}\n"
                          "No. order reds due to stab lim detection   {}\n"
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

  // ------------ Direct solvers interface

  Direct::Direct()
    : ICVODEStrategy()
    , matrix{ nullptr }
    {}
  
  Direct::~Direct() {
    if( matrix ) SUNMatDestroy(matrix);
  }
  
  int Direct::BindUserFunctions( ) {
    return ( ( client->opt_udf_set.jacobian ) 
             ? CVodeSetJacFn( memory, Client::JacobianCallback )
             : 0 );
  }

  // ------------ Dense matrix direct strategy
  
  Dense::Dense()
    : ICVODEStrategy()
    , Direct()
  {}
  
  Dense::~Dense() {}
  
  int Dense::SetLinearSolver() {
    // create dense SUNMatrix for use in linear solves
    matrix = SUNDenseMatrix( client->n_states, client->n_states );
    // create the linear solver
    linear_solver = SUNLinSol_Dense( client->state, matrix );
    return CVodeSetLinearSolver( memory
                                 , linear_solver
                                 , matrix );
  }
  
  // ------------ Band matrix direct strategy
  
  Band::Band()
    : ICVODEStrategy()
    , Direct()
    , bandwidth()
  {}
    
  Band::~Band() {}
  
  void Band::SetMatrixBandwidth( Types::Bandwidth const bandwidth_ ) {
    bandwidth = bandwidth_;
  }
  
  void Band::SetMatrixBandwidth( sunindextype const upper
                                 , sunindextype const lower ) {
    bandwidth.upper = upper;
    bandwidth.lower = lower;
  }
  
  int Band::SetLinearSolver() {
    // create band SUNMatrix for use in linear solves
    assert( bandwidth.upper >= 0 );
    assert( bandwidth.lower >= 0 );
    matrix = SUNBandMatrix( client->n_states
                            , bandwidth.upper
                            , bandwidth.lower );
    // create the linear solver
    linear_solver = SUNLinSol_Band( client->state, matrix );
    return CVodeSetLinearSolver( memory
                                 , linear_solver
                                 , matrix );
  }
  
    
  // ------------ Iterative solvers interface
  
  Iterative::Iterative()
    : ICVODEStrategy()
    , options()
  { }
  
  Iterative::~Iterative() {}
  
  void Iterative::SetLinearSolverOptions( Types::Preconditioner const preconditioner
                                          , Types::GramSchmidt const gram_schmidt
                                          , int const n_krylov_basis_vectors
                                          , int const max_restarts ) {
    options = { preconditioner
                , gram_schmidt
                , n_krylov_basis_vectors
                , max_restarts };
  }
             
  int Iterative::BindUserFunctions() {
    // set the Jacobian functions
    int ret_code;
    if( client->opt_udf_set.jacobian_times_vector ) {
      CVLsJacTimesSetupFn jacobian_times_vector_setup_func =
        ( client->opt_udf_set.jacobian_times_vector_setup )
        ? Client::JacobianTimesVectorSetupCallback
        : nullptr;
      ret_code = CVodeSetJacTimes( ICVODEStrategy::memory
                                   , jacobian_times_vector_setup_func
                                   , Client::JacobianTimesVectorCallback );
      if( ret_code != 0 ) return ret_code;
    }
    
    // Set the preconditioner functions
    ret_code = CVodeSetPreconditioner( memory
                                       , Client::PreconditionerSetupCallback
                                       , Client::PreconditionerSolveCallback );
    if( ret_code != 0 ) return ret_code;
    
    return ret_code;
  }

  // ------------ SPGMR iterative strategy
  
  SPGMR::SPGMR()
    : ICVODEStrategy()
    , Iterative()
  { }
  
  SPGMR::~SPGMR() { }
  
  int SPGMR::SetLinearSolver() {
    // create the linear solver
    linear_solver = SUNLinSol_SPGMR( client->state
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

  // ------------ SPFGMR iterative strategy

  SPFGMR::SPFGMR()
    : ICVODEStrategy()
    , Iterative()
  { }
  
  SPFGMR::~SPFGMR() { }
     
  int SPFGMR::SetLinearSolver() {
    // create the linear solver
    linear_solver = SUNLinSol_SPFGMR( client->state
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

  // ------------ SPBCGS iterative strategy
     
  SPBCGS::SPBCGS()
    : ICVODEStrategy()
    , Iterative()
  { }
  
  SPBCGS::~SPBCGS() { }
    
  int SPBCGS::SetLinearSolver() {
    // create the linear solver
    linear_solver = SUNLinSol_SPBCGS( client->state
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

  // ------------ SPTFQMR iterative strategy
  
  SPTFQMR::SPTFQMR()
    : ICVODEStrategy()
    , Iterative()
  { }
  
  SPTFQMR::~SPTFQMR() { }
  
  int SPTFQMR::SetLinearSolver() {
    // create the linear solver
    linear_solver = SUNLinSol_SPTFQMR( client->state
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
  
  
  
} // namespace SUNDIALS
