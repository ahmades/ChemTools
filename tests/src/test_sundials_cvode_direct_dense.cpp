#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>

#include <cmath>

#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/color.h>

#include <catch2/catch.hpp>

#include "test_config.hpp"
#include "sundials_cvode.hpp"

struct Time {
  realtype start;
  realtype stop;
  Time( realtype start_
        ,  realtype stop_ )
    : start( start_ )
    , stop( stop_ )
  {}
};

class CVRobertsDNS : public SUNDIALS::Client {
public:

  CVRobertsDNS( std::vector<realtype> const& state
                , Time const & time
                , SUNDIALS::Types::IntegrationTolerance const& tolerance )
    : m_state( state )
    , m_time( time )
    , m_tolerance( tolerance )
  {
    // enable optional user functions
    Client::opt_udf_set.jacobian = true;
  }

  void Initialise() {
    // set initial conditions
    SetState( m_state );
    // set tolerances: optional
    SetTolerances( m_tolerance.relative, m_tolerance.absolute );
    // integration time
    SetIntegrationTime( m_time.start, m_time.stop );
  }

  int RightHandSide( realtype const /*time*/
                     , realtype* const state
                     , realtype* rhs ) override {
    
    rhs[0] = -0.04 * state[0] + 1.0e+4* state[1] * state[2];
    rhs[2] = 3.0e+7 * state[1] * state[1];
    rhs[1] = -rhs[0] - rhs[2];
    
    return 0;
  }

  int Jacobian( realtype const /*time*/
                , realtype* const state
                , realtype* const /*rhs*/
                , realtype** jacobian_mat ) override {

    jacobian_mat[0][0] = -0.04;
    jacobian_mat[1][0] = 1.0e+4 * state[2];
    jacobian_mat[2][0] = 1.0e+4 * state[1];
    
    jacobian_mat[0][1] = 0.04; 
    jacobian_mat[1][1] = -1.0e+4 * state[2] - 6.0e+7 * state[1];
    jacobian_mat[2][1] = -1.0e+4 * state[1];
    
    jacobian_mat[0][2] = 0.0;
    jacobian_mat[1][2] = 6.0e+7 * state[1];
    jacobian_mat[2][2] = 0.0;
    
    return 0;
  }

private:
  
  std::vector<realtype> m_state;                     // input: initial state
  Time m_time;                                       // input: time info
  SUNDIALS::Types::IntegrationTolerance m_tolerance; // input: tolerances
  
};


TEST_CASE( "3-species reaction PDE system can be solved", "[direct dense solver]" ) {

  std::vector<realtype> state_init{ 1.0, 0.0, 0.0 };
  
  Time const time( 0.0, std::numeric_limits<double>::infinity() );
  
  realtype const rel_tol = 1.0e-4; // scalar relative tolerance
  std::vector<realtype> const abs_tol{ 1.0e-8, 1.0e-14, 1.0e-6}; // vector absolute tolerance
  SUNDIALS::Types::IntegrationTolerance const tolerance( rel_tol, abs_tol );
  
  CVRobertsDNS client( state_init
                       , time
                       , tolerance );

  // Initialise the client
  client.Initialise();  

  // select strategy
  SUNDIALS::Dense strategy;
  // set strategy
  SUNDIALS::CVODE cvode( &strategy );
  
  // set linear multistep method
  cvode.SetLinearMultiStepMethod( SUNDIALS::Types::LinearMultisptepMethod::BDF );
      
  // set client data
  cvode.SetClientData( &client );

 // set time step control options
  cvode.SetInitStepSize( 0.0 ); // 0 sets default
  cvode.SetMinStepSize( 0.0 );  // 0 sets default
  cvode.SetMaxStepSize( 0.0 );  // 0 sets default
  cvode.SetMaxNumSteps( 0 );    // 0 sets default, 500. 600 > 500 => should not influence the test

  // set solver control options
  cvode.SetMaxOrder( 5 );                        // same value as target, otherwise test will fail
  cvode.SetMaxWarnMessages( 15 );                // target uses default, 10, does not influence the test
  cvode.SetStabilityLimitDetection( true );
  cvode.SetMaxErrorTestFailures( 7 );            // target uses default, 7, but test should pass
  cvode.SetMaxNonlinearIterations( 3 );          // target uses default, 3, but test should pass
  cvode.SetMaxConvergenceFailures(15);           // target uses default, 10, but test should pass
  cvode.SetNonlinConvergenceCoefficient( 0.1 );  // same value as target, otherwise test will fail
  
  // initialise
  REQUIRE( 0 == cvode.Initialise() );

  // integrate
  {
    realtype time = 0.4;
    realtype const time_step_mult = 10.0;
    size_t const n_time_steps = 12;
    size_t time_step_count = 0;
    std::vector<realtype> test_state;         // vector for storage of temporal solutions
    
    while( time_step_count < n_time_steps ) {
      // integrate
      REQUIRE( CV_SUCCESS == cvode.Integrate( time ) );

      // append the vector of states
      test_state.insert( test_state.end()
                         , client.State()
                         , client.State() + client.NStates() );
      
      // log if verbose
      if( TestConfig::verbose ) {
        fmt::print( "At t = {:e} s:  y[0] = {:e}  y[1] = {:e}  y[2] = {:e}\n"
                    , time
                    , client.State( 0 )
                    , client.State( 1 )
                    , client.State( 2 ) );
      }
      
      time_step_count++;
      time *= time_step_mult;      
    }
    
    if( TestConfig::verbose ) {
      std::string const title = fmt::format( fmt::emphasis::underline
                                             | fmt::emphasis::bold
                                             , "Final statistics:" );    
      fmt::print( "{}\n{}\n"
                  , title
                  , cvode.PrintSolverStatistics() );
    }
    
    // read and store expected state from test data dir
    std::vector<realtype> expected_state;
    {
      std::ifstream file;
      realtype val;
      file.open( "../tests/data/sundials_test_cvRoberts_dns.dat", std::ifstream::in );
      if( file.is_open() ) {
        while ( file >> val )  {
          expected_state.push_back( val );
        }
        file.close();
      }
    }
    
    // computed and expected state vectors must have the same size
    REQUIRE( test_state.size() == expected_state.size());
    
    // computed and expected state vectors must be equal within a tiny tolerance
    REQUIRE_THAT( test_state
                  , Catch::Approx( expected_state ).epsilon( TestConfig::tolerance ) );
    
  }
  
}
