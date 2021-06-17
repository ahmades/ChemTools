#ifndef SUNDIALS_CVODE_INTERFACE_HPP
#define SUNDIALS_CVODE_INTERFACE_HPP

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

#include <boost/optional.hpp>

#include "sundials/cvode/types.hpp"
#include "sundials/cvode/client.hpp"

namespace SUNDIALS {

  namespace CVODE {
  
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
    class IStrategy {

    protected:
    
      void* memory;
      Client* client;
      SUNLinearSolver linear_solver;
    
    private:
    
      Types::LinearMultisptepMethod linear_multi_step_meth;
      Types::TimeStepcontrol time_step_ctrl;
      Types::SolverControl solver_ctrl;
      bool is_initialised;
      bool must_be_reinitialised;
    
      virtual int SetLinearSolver() = 0;
      virtual int BindUserFunctions() = 0;
    
    public:

      IStrategy();
      virtual ~IStrategy();
      void SetLinearMultiStepMethod( Types::LinearMultisptepMethod const linear_multi_step_meth_ );
      void SetClientData( void* const user_data_ );
      boost::optional<int> SetInitStepSize( realtype const init_step );
      boost::optional<int> SetMinStepSize( realtype const min_step );
      boost::optional<int> SetMaxStepSize( realtype const max_step );
      boost::optional<int> SetMaxNumSteps( long int const n_steps );
      void SetMaxOrder( int const max_order );
      boost::optional<int> SetMaxWarnMessages( int const max_warn_msgs );
      boost::optional<int> SetStabilityLimitDetection( bool const stab_lib_det_active );
      boost::optional<int> SetMaxErrorTestFailures( int const max_err_test_fails );
      boost::optional<int> SetMaxNonlinearIterations( int const max_nonlin_iters );
      boost::optional<int> SetMaxConvergenceFailures( int const max_conv_fails );
      boost::optional<int> SetNonlinConvergenceCoefficient( realtype const nonlin_conv_coef );  
      int Initialise();
      int ReInitialise( realtype const time_init,
                        std::vector<realtype> const& state );
      int ReInitialise( realtype const time_init,
                        realtype const * const state );  
      int Integrate( realtype const time_target );
      realtype* Solution() const;
      realtype Solution( sunindextype const i ) const;
      int MainSolverStatistics( Types::MainSolverStatistics& stats ) const;
      int LinearSolverStatistics( Types::LinearSolverStatistics& stats ) const;
      std::string PrintSolverStatistics() const;
    };

    // ------------ Direct solvers interface
    class Direct: public virtual IStrategy {
    public:
      Direct();
      ~Direct();        
    
    protected:
      SUNMatrix matrix;

    private:
      int BindUserFunctions() override;
    
    };

    // ------------ Dense matrix direct strategy
    class Dense : public virtual IStrategy
                , public Direct {
    public:
      Dense();    
      ~Dense();
    
    private:
      int SetLinearSolver() override;
    };

    // ------------ Band matrix direct strategy
    class Band : public virtual IStrategy
               , public Direct {
    public:
      Band();
      ~Band();
      void SetMatrixBandwidth( Types::Bandwidth const bandwidth_ );
      void SetMatrixBandwidth( sunindextype const upper
                               , sunindextype const lower );

    private:
      Types::Bandwidth bandwidth;
      int SetLinearSolver() override;    
    };
 
    // ------------ Iterative solvers interface
    class Iterative: public virtual IStrategy {
    public:
      Iterative();    
      ~Iterative();
      void SetLinearSolverOptions( Types::Preconditioner const preconditioner
                                   , Types::GramSchmidt const gram_schmidt
                                   , int const n_krylov_basis_vectors
                                   , int const max_restarts );
  
    private:
      int BindUserFunctions() override;

    protected:
      Types::IterativeLinarSolverOptions options;
    };

    // ------------ SPGMR iterative strategy
    class SPGMR : public virtual IStrategy
                , public Iterative {
    public:
      SPGMR();
      ~SPGMR();
    
    private:
      int SetLinearSolver() override;
    };

    // ------------ SPFGMR iterative strategy
    class SPFGMR : public virtual IStrategy
                 , public Iterative {
    public:
      SPFGMR();
      ~SPFGMR();
    
    private:
      int SetLinearSolver() override;
    };

    // ------------ SPBCGS iterative strategy
    class SPBCGS : public virtual IStrategy
                 , public Iterative {
    public:
      SPBCGS();
      ~SPBCGS();
    
    private:
      int SetLinearSolver() override;
    };
  
    // ------------ SPTFQMR iterative strategy
    class SPTFQMR : public virtual IStrategy
                  , public Iterative {
    public:
      SPTFQMR();    
      ~SPTFQMR();
    
    private:
      int SetLinearSolver() override;
    };
  

    // CVODE client - context
    class Solver {
    
    private:
    
      IStrategy* strategy;

    public:

      explicit Solver()
        : strategy( nullptr )
      {}
    
      explicit Solver( IStrategy* strategy_ )
        : strategy( strategy_ )
      {}
    
      void SetStrategy( IStrategy* strategy_ ) {
        strategy = strategy_;
      }

      void SetLinearMultiStepMethod( Types::LinearMultisptepMethod
                                     const linear_multi_step_meth ) {
        strategy->SetLinearMultiStepMethod( linear_multi_step_meth );
      }
    
      void SetClientData( void* const user_data ) {
        strategy->SetClientData( user_data );
      }

      boost::optional<int> SetInitStepSize( realtype const init_step ) const {
        return strategy->SetInitStepSize( init_step );
      }
    
      boost::optional<int> SetMinStepSize( realtype const min_step ) const {
        return strategy->SetMinStepSize( min_step );
      }
    
      boost::optional<int> SetMaxStepSize( realtype const max_step ) const {
        return strategy->SetMaxStepSize( max_step );
      }
    
      boost::optional<int> SetMaxNumSteps( long int const n_steps ) const {
        return strategy->SetMaxNumSteps( n_steps );
      }

      void SetMaxOrder( int const max_order ) const {
        return strategy->SetMaxOrder( max_order );
      }
    
      boost::optional<int> SetMaxWarnMessages( int const max_warn_msgs ) const {
        return strategy->SetMaxWarnMessages( max_warn_msgs );
      }
    
      boost::optional<int> SetStabilityLimitDetection( bool const stab_lib_det_active ) const {
        return strategy->SetStabilityLimitDetection( stab_lib_det_active );
      }
    
      boost::optional<int> SetMaxErrorTestFailures( int const max_err_test_fails ) const {
        return strategy->SetMaxErrorTestFailures( max_err_test_fails );
      }
    
      boost::optional<int> SetMaxNonlinearIterations( int const max_nonlin_iters ) const {
        return strategy->SetMaxNonlinearIterations( max_nonlin_iters );
      }
    
      boost::optional<int> SetMaxConvergenceFailures( int const max_conv_fails ) const {
        return strategy->SetMaxConvergenceFailures( max_conv_fails );
      }
    
      boost::optional<int> SetNonlinConvergenceCoefficient( realtype const
                                                            nonlin_conv_coef ) const {
        return strategy->SetNonlinConvergenceCoefficient( nonlin_conv_coef );
      }
    
      int Initialise() const {
        return strategy->Initialise();
      }
    
      int Integrate( realtype const time_target ) const {
        return strategy->Integrate( time_target );
      }

      realtype* Solution() const {
        return strategy->Solution( );
      }
    
      int MainSolverStatistics( Types::MainSolverStatistics& stats ) const {
        return strategy->MainSolverStatistics( stats );
      }

      int LinearSolverStatistics( Types::LinearSolverStatistics& stats ) const {
        return strategy->LinearSolverStatistics( stats );
      }
    
      std::string PrintSolverStatistics() const {
        return strategy->PrintSolverStatistics();
      }
    };

  } // namespace CVODE

} // namespace SUNDIALS


#endif // SUNDIALS_CVODE_INTERFACE_HPP
