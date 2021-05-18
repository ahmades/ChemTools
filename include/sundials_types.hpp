#ifndef SUNDIALS_TYPES_HPP
#define SUNDIALS_TYPES_HPP

#include <cvode/cvode.h>               // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>    // access to serial N_Vector  
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix 
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sundials/sundials_types.h>   // defs. of realtype, sunindextype

namespace SUNDIALS {

  namespace Types {

    // cvode types

    enum class LinearMultisptepMethod {
        invalid = 0
      , Adams   = CV_ADAMS // 1
      , BDF     = CV_BDF   // 2
    };

    enum class Preconditioner {
        none  = PREC_NONE
      , left  = PREC_LEFT
      , right = PREC_RIGHT
      , both  = PREC_BOTH
    };

    enum class GramSchmidt {
        modified  = MODIFIED_GS
      , classical = CLASSICAL_GS
    };

    enum class Constraint : int {
        none              =  0
      , positive          =  1
      , negative          = -1
      , strictly_positive =  2
      , strictly_negative = -2
    };
    
    class IntegrationTime {
    public:
      realtype start;
      realtype stop;
      realtype step;
      realtype reached;
      
      IntegrationTime( )
        : start{ 0.0 }
        , stop{ std::numeric_limits<realtype>::infinity( ) }
        , step{ 0.0 }
        , reached{ 0.0 }
      { }
    };

    class IntegrationTolerance {
    public:
      realtype relative;
      realtype absolute;
      
      IntegrationTolerance( )
        : relative{ 1.0e-6 }
        , absolute{ 1.0e-6 }
      { }
    };

    class CommonUserFunctions {
    public:
      CVRhsFn right_hand_side;
      CVProjFn projection;
      CVErrHandlerFn error_handler;
      
      CommonUserFunctions( )
        : right_hand_side{ nullptr }
        , projection{ nullptr }
        , error_handler{ nullptr }
      { }
    };

    class DirectUserFunctions {
    public:
      CVLsJacFn jacobian;
      
      DirectUserFunctions( )
        : jacobian{ nullptr }
      { }
    };
    
    class IterativeUserFunctions {
    public:
      CVLsJacTimesSetupFn jacobian_setup;
      CVLsJacTimesVecFn jacobian_times_vector;
      CVLsPrecSetupFn preconditioner_setup;
      CVLsPrecSolveFn preconditioner_solve;
      
      IterativeUserFunctions( )
        : jacobian_setup{ nullptr }
        , jacobian_times_vector{ nullptr }
        , preconditioner_setup{ nullptr }
        , preconditioner_solve{ nullptr }
      { }
    };

    class IterativeLinarSolverOptions {
    public:
      Preconditioner preconditioner;
      GramSchmidt gram_schmidt;
      int n_krylov_basis_vectors;
      int max_restarts;
      
      IterativeLinarSolverOptions()
        : preconditioner{ Preconditioner::none }
        , gram_schmidt{ GramSchmidt::modified }
        , n_krylov_basis_vectors{ 5 }
        , max_restarts{ 0 }
      {}

      IterativeLinarSolverOptions( Preconditioner const preconditioner_
                                   , GramSchmidt const gram_schmidt_
                                   , int const n_krylov_basis_vectors_
                                   , int const max_restarts_ )
        : preconditioner{ preconditioner_ }
        , gram_schmidt{ gram_schmidt_ }
        , n_krylov_basis_vectors{ n_krylov_basis_vectors_ }
        , max_restarts{ max_restarts_ }
      {}
    };

    class Bandwidth {
    public:
      sunindextype upper;
      sunindextype lower;
    };
   
    class MainSolverStatistics {
    public:
      long int n_internal_steps;
      long int n_rhs_evals;
      long int n_linear_solver_setups;
      long int n_error_test_fails;
      long int n_nonlinear_solver_iters;
      long int n_nonlinear_solver_conv_fails;
      int last_order_used;
      int current_order;
      realtype last_internal_step_size ;
      realtype next_internal_step_size ;
      realtype first_internal_step_size;
      long int n_stability_order_reductions;
      MainSolverStatistics()
        : n_internal_steps{ 0 }
        , n_rhs_evals{ 0 }
        , n_linear_solver_setups{ 0 }
        , n_error_test_fails{ 0 }
        , n_nonlinear_solver_iters{ 0 }
        , n_nonlinear_solver_conv_fails{ 0 }
        , last_order_used{ 0 }
        , current_order{ 0 }
        , last_internal_step_size { 0.0 }
        , next_internal_step_size { 0.0 }
        , first_internal_step_size{ 0.0 }
        , n_stability_order_reductions{ 0 }
      {}
    };

    class LinearSolverStatistics {
    public:
      long int n_jac_evals;
      long int n_rhs_jac_evals;
      long int n_linear_iterations;
      long int n_linear_conv_fails;
      long int n_prec_evals;
      long int n_prec_solves;
      long int n_jac_vec_setup_evals;
      long int n_jac_vec_prod_evals;
      LinearSolverStatistics()
        : n_jac_evals{ 0 }
        , n_rhs_jac_evals{ 0 }
        , n_linear_iterations{ 0 }
        , n_linear_conv_fails{ 0 }
        , n_prec_evals{ 0 }
        , n_prec_solves{ 0 }
        , n_jac_vec_setup_evals{ 0 }
        , n_jac_vec_prod_evals{ 0 }
      {}
    };
    
  } // namespace Types
} // namespace SUNDIALS

#endif //SUNDIALS_TYPES_HPP
