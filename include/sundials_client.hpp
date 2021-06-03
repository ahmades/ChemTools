#ifndef SUNDIALS_CLIENT_HPP
#define SUNDIALS_CLIENT_HPP

#include <vector>
#include <algorithm>
#include <stdexcept>

#include <nvector/nvector_serial.h>

#include "sundials_types.hpp"

namespace SUNDIALS {

  /* 
     CVODE client base class:
      - All classes that use CVODE must inherit from Client.
      - All derived classes must override a set of virtual methods, depending on the strategy.
      - members of Sub-class OptUDFSet must be set to true accordingly.
        TODO: this is bad, must find a better way.
   */
  
  class Client {
   
  public:

    size_t n_states;     // number of states
    N_Vector state;      // state vector
    N_Vector constraint; // constraints vector
    bool has_constraints;
    Types::IntegrationTime time;
    Types::IntegrationTolerance tolerance;

    // switches for optional user-defined functions
    class OptUDFSet {
    public:
      bool jacobian;
      bool jacobian_times_vector_setup;
      bool jacobian_times_vector;
      bool projection;
      bool error_handler;
      OptUDFSet()
        : jacobian{ false }
        , jacobian_times_vector_setup{ false }
        , jacobian_times_vector{ false }
        , projection{ false }
        , error_handler{ false }
      { }
    } opt_udf_set;

    // constructor
    Client()
      : n_states{ 0 }
      , state{ nullptr }
      , constraint{ nullptr }
      , has_constraints{ false }
      , time()
      , tolerance()
      , opt_udf_set( )
    {}

    // destructor
    virtual ~Client() {
      if( state ) N_VDestroy_Serial( state );
      if( constraint ) N_VDestroy_Serial( constraint );
    }

    // allocates and sets state and optional constraints vectors
    void SetState( std::vector<realtype> const& state_vec
                   , std::vector<Types::Constraint> const& constraint_vec
                   = std::vector<Types::Constraint>() ) {
      n_states = state_vec.size();
      assert( n_states > 0 );
      state = N_VNew_Serial( static_cast<sunindextype>( n_states ) );
      realtype * const state_ptr = N_VGetArrayPointer( state );
      std::copy( state_vec.begin()
                 , state_vec.end()
                 , state_ptr );
      if( !constraint_vec.empty() ) {
        assert( n_states == constraint_vec.size() );
        has_constraints = true;
        constraint = N_VNew_Serial( static_cast<sunindextype>( n_states ) );
        std::vector<realtype> tmp;
        std::transform( constraint_vec.begin()
                        , constraint_vec.end()
                        , std::back_inserter( tmp )
                        , []( Types::Constraint const con )
                        { return static_cast<double>( con ); } );
        realtype* const constraints_ptr = N_VGetArrayPointer( constraint );
        std::copy( tmp.begin()
                   , tmp.end()
                   , constraints_ptr );
      }
    }

    size_t NStates() const {
      return n_states;
    };

    realtype* State() const {
      assert( state );
      return N_VGetArrayPointer( state );
    }

    realtype State( sunindextype const i ) const {
      assert( state );
      return N_VGetArrayPointer( state )[i];
    }

    void SetTolerances( realtype const relative_tolerance
                        , realtype const absolute_tolerance ) {
      tolerance.relative = relative_tolerance;
      tolerance.absolute.front() = absolute_tolerance;
    }

    void SetTolerances( realtype const relative_tolerance
                        , std::vector<realtype> const& absolute_tolerance ) {
      tolerance.relative = relative_tolerance;
      tolerance.absolute = absolute_tolerance;
    }
  
    void SetIntegrationTime( realtype const time_start
                             , realtype const time_stop ) {
      assert( time_start >= 0.0 );
      assert( time_stop > time_start );
      time.start = time_start;
      time.stop = time_stop;
    }

    // Callback functions

    static int RightHandSideCallback( realtype const time
                                      , N_Vector const state_nvec
                                      , N_Vector rhs_nvec
                                      , void* const client_ ) {
      Client* const client = static_cast<Client* const>( client_ );
      return client->RightHandSide( time
                                    , N_VGetArrayPointer( state_nvec )
                                    , N_VGetArrayPointer( rhs_nvec ) );
    }
    
    static int JacobianCallback( realtype const time
                                 , N_Vector const state_nvec
                                 , N_Vector const rhs_nvec
                                 , SUNMatrix jacobian_mat
                                 , void* const client_
                                 , N_Vector /*tmp_1*/
                                 , N_Vector /*tmp_2*/
                                 , N_Vector /*tmp_3*/ ) {
      Client* const client = static_cast<Client* const>( client_ );
      return client->Jacobian( time
                               , N_VGetArrayPointer( state_nvec )
                               , N_VGetArrayPointer( rhs_nvec )
                               , (SUNMATRIX_DENSE == SUNMatGetID( jacobian_mat) )
                               ? SUNDenseMatrix_Cols( jacobian_mat )
                               : SUNBandMatrix_Cols( jacobian_mat ) );
    }
    
    static int JacobianTimesVectorSetupCallback( realtype const time
                                                 , N_Vector const state_nvec
                                                 , N_Vector const rhs_nvec
                                                 , void* const client_) {
      Client* const client = static_cast<Client* const>( client_ );
      return client->JacobianTimesVectorSetup( time
                                               , N_VGetArrayPointer( state_nvec )
                                               , N_VGetArrayPointer( rhs_nvec ) );
    }
    
    
    static int JacobianTimesVectorCallback( N_Vector const vector_nvec
                                            , N_Vector jacobian_vector_nvec
                                            , realtype const time
                                            , N_Vector const state_nvec
                                            , N_Vector const rhs_nvec
                                            , void* const client_
                                            , N_Vector /*tmp*/ ) {
      Client* const client = static_cast<Client* const>( client_ );
      return client->JacobianTimesVector( N_VGetArrayPointer( vector_nvec )
                                          , N_VGetArrayPointer( jacobian_vector_nvec )
                                          , time
                                          , N_VGetArrayPointer( state_nvec )
                                          , N_VGetArrayPointer( rhs_nvec ) );
    }
    
    static int PreconditionerSetupCallback( realtype const time
                                            , N_Vector const state_nvec
                                            , N_Vector const rhs_nvec
                                            , booleantype const jac_ok
                                            , booleantype* jac_cur_ptr
                                            , realtype const gamma
                                            , void* const client_ ) {
      Client* const client = static_cast<Client* const>( client_ );
      return client->PreconditionerSetup( time
                                          , N_VGetArrayPointer( state_nvec )
                                          , N_VGetArrayPointer( rhs_nvec )
                                          , jac_ok
                                          , jac_cur_ptr
                                          , gamma );
    }
    
    static int PreconditionerSolveCallback( realtype const time
                                            , N_Vector const state_nvec
                                            , N_Vector const rhs_nvec
                                            , N_Vector const r_nvec
                                            , N_Vector z_nvec
                                            , realtype const gamma
                                            , realtype const delta
                                            , int const lr
                                            , void* const client_ ) {
      Client* const client = static_cast<Client* const>( client_ );
      return client->PreconditionerSolve( time
                                          , N_VGetArrayPointer( state_nvec )
                                          , N_VGetArrayPointer( rhs_nvec )
                                          , N_VGetArrayPointer( r_nvec )
                                          , N_VGetArrayPointer( z_nvec )
                                          , gamma
                                          , delta
                                          , lr );
    }
    
    static void ErrorHandlerCallback( int const error_code
                                      , const char* module
                                      , const char* function
                                      , char* message
                                      , void* const client_ ) {
      Client* const client = static_cast<Client* const>( client_ );
      client->ErrorHandler( error_code
                            , module
                            , function
                            , message );
    }
    
    static int ProjectionCallback( realtype const time
                                   , N_Vector const state_nvec
                                   , N_Vector const correction_nvec
                                   , realtype const projection_tolerane
                                   , N_Vector projection_error_nvec
                                   , void* const client_ ) {
      Client* const client = static_cast<Client* const>( client_ );
      return client->Projection( time
                                 , N_VGetArrayPointer( state_nvec )
                                 , N_VGetArrayPointer( correction_nvec )
                                 , projection_tolerane
                                 , N_VGetArrayPointer( projection_error_nvec ) );
    }

  protected:
    
    // Right-hand-side function
    virtual int RightHandSide( realtype const time
                               , realtype* const state
                               , realtype* rhs ) = 0;
    
    // Optional Jacobian function for direct solvers   
    virtual int Jacobian( realtype const /*time*/
                          , realtype* const /*state*/
                          , realtype* const /*rhs*/
                          , realtype**/* jacobian_mat*/ ) {
      throw std::runtime_error( "SUNDIALS::Client: a direct solver was selected "
                                "but the user-defined function Jacobian was not supplied." );
      return -1;
    }

    // Optional Jacobian times vector setup function for iterative solvers
    virtual int JacobianTimesVectorSetup( realtype const /*time*/
                                          , realtype* const /*state*/
                                          , realtype* const /*rhs*/ ) {
      throw std::runtime_error( "SUNDIALS::Client: an iterative solver was selected "
                                "but the optional user-defined function JacobianTimesVectorSetup "
                                "was not supplied." );
      return -1;
    }

    // Optional Jacobian times vector function for iterative solvers
    virtual int JacobianTimesVector( realtype* const/* vector*/
                                     , realtype* /*jacobian_vector*/
                                     , realtype const /*time*/
                                     , realtype* const /*state*/
                                     , realtype* const/* rhs*/ ) {
      throw std::runtime_error( "SUNDIALS::Client: an iterative solver was selected "
                                "but the optional user-defined function JacobianTimesVector "
                                "was not supplied." );
      return -1;
    }

    // Preconditioner setup function for iterative solvers
    virtual int PreconditionerSetup( realtype const /*time*/
                                     , realtype* const /*state*/
                                     , realtype* const/* rhs*/
                                     , booleantype const /*jac_ok*/
                                     , booleantype* /*jac_cur_ptr*/
                                     , realtype const /*gamma*/ ) {
      throw std::runtime_error( "SUNDIALS::Client: an iterative solver was selected "
                                "but the user-defined function PreconditionerSetup was not supplied." );
      return -1;
    }

    // Preconditioner solve function for iterative solvers
    virtual int PreconditionerSolve( realtype const /*time*/
                                     , realtype* const /*state*/
                                     , realtype* const /*rhs*/
                                     , realtype* const /*r*/
                                     , realtype* /*z*/
                                     , realtype const /*gamma*/
                                     , realtype const /*delta*/
                                     , int const /*lr*/ )  {
      throw std::runtime_error( "SUNDIALS::Client: an iterative solver was selected "
                                "but the user-defined function PreconditionerSolve was not supplied." );
      return -1;
    }

    // Optional projection function
    virtual int Projection( realtype const /*time*/
                            , realtype* const /*state*/
                            , realtype* const /*correction*/
                            , realtype const /*projection_tolerane*/
                            , realtype* /*projection_error*/ )  {
      throw std::runtime_error( "SUNDIALS::Client: a constrained solution was requested "
                                "but the user-defined function Projection was not supplied." );
      return -1;
    }

    // Optional error handler function
    virtual void ErrorHandler( int const /*error_code*/
                               , const char* /*module*/
                               , const char* /*function*/
                               , char* /*message*/ ) {
      throw std::runtime_error( "SUNDIALS::Client: a custom error handler was requested "
                                "but the user-defined function ErrorHandler was not supplied." );
      return;
    }
  };
  
} // namespace SUNDIALS

#endif // SUNDIALS_CLIENT_HPP
