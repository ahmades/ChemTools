#include <vector>
#include <limits>
#include <fstream>

#include <cmath>

#include <fmt/format.h>
#include <catch2/catch.hpp>

#include "test_config.hpp"
#include "sundials_cvode.hpp"


/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 *                
 * Apart from the test, modied by yours truly to have a more c++y feel:
 *  - got rid of most macros
 *  - encapsulated the whole thing in a single class
 *  - improved formatting and readability
 *  - applied some arithmetic optimisations
 *  - enforced const correctness
 *
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(y) = Kv0*exp(y/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * 10 x 10 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODE, with the BDF/GMRES
 * method (i.e. using the SUNLinSol_SPGMR linear solver) and the
 * block-diagonal part of the Newton matrix as a left
 * preconditioner. A copy of the block-diagonal part of the
 * Jacobian is saved and conditionally reused within the Precond
 * routine.
 * -----------------------------------------------------------------*/

class CVDiurnalKry;
#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*CVDiurnalKry::NS + (k)*CVDiurnalKry::NSNX])

class CVDiurnalKry {
public:

  static realtype constexpr X_MIN     = 0.0;                     // grid boundaries in x: min
  static realtype constexpr X_MAX     = 20.0;                    // grid boundaries in x: max
  static realtype constexpr Y_MIN     = 30.0;                    // grid boundaries in y: min
  static realtype constexpr Y_MAX     = 50.0;                    // grid boundaries in y: max
  static realtype constexpr X_MID     = ( X_MIN + X_MAX ) / 2.0; // grid midpoints in x
  static realtype constexpr Y_MID     = ( Y_MIN + Y_MAX ) / 2.0; // grid midpoints in y
  
  static int      constexpr NS        = 2;                       // number of species
  static int      constexpr NX        = 10;                      // number of x mesh points 
  static int      constexpr NY        = 10;                      // number of y mesh points 
  static int      constexpr NSNX      = NS*NX;                   // NSNX = NS*NX 
  static int      constexpr NEQ       = NS*NX*NY;                // total number of equations 
  
  static realtype constexpr KH        = 4.0e-6;                  // horizontal diffusivity Kh 
  static realtype constexpr VEL       = 0.001;                   // advection velocity V      
  static realtype constexpr KV        = 1.0e-8;                  // coefficient in Kv;      
  static realtype constexpr COEF_Q_1  = 1.63e-16;                // coefficient q1
  static realtype constexpr COEF_Q_2  = 4.66e-16;                // coefficient q2
  static realtype constexpr COEF_C_3  = 3.7e+16;                 // coefficient c3   
  static realtype constexpr COEF_A_3  = 22.62;                   // coefficient in expression for q3(t) 
  static realtype constexpr COEF_A_4  = 7.601;                   // coefficient in expression for q4(t)
  static realtype constexpr SCALE_C_1 = 1.0e6;                   // coefficients in initial profiles    
  static realtype constexpr SCALE_C_2 = 1.0e+12;

  static realtype constexpr T0       = 0.0;                      // initial time
  static int      constexpr NOUT     = 12;                       // number of output times
  static realtype constexpr TWO_HR   = 7200.0;                   // number of seconds in two hours 
  static realtype constexpr HALF_DAY = 4.32e+4;                  // number of seconds in a half day
  static realtype constexpr PI       = 3.1415926535898;          // pi (not using M_PI from cmath to match sundial's results
  
  static realtype constexpr REL_TOL  = 1.0e-5;                   // scalar relative tolerance
  static realtype constexpr ABS_TOL  = 100.0 * REL_TOL;          // scalar absolute tolerance

  static realtype* Access( realtype* data, int i, int j, int k ) {
    return &data[i-1 + j*NS + k*NSNX];
  }

  // static realtype& Access( std::vector<realtype>& vdata, int i, int j, int k ) {
  //   return vdata[i-1 + j*NS + k*NSNX];
  // }
   
  class Data {
  public:
    realtype q_4;
    realtype om;
    realtype dx;
    realtype dy;
    realtype hdco;
    realtype haco;
    realtype vdco;
    realtype **P[NX][NY];
    realtype **Jbd[NX][NY];
    sunindextype *pivot[NX][NY];

    Data()
      : q_4{ 0.0 }
      , om{ PI / HALF_DAY }
      , dx{ ( X_MAX - X_MIN ) / ( NX - 1) }
      , dy{ ( Y_MAX - Y_MIN ) / ( NY - 1 ) }
      , hdco{ KH/std::pow( dx, 2 ) }
      , haco{ VEL / ( 2.0 * dx ) }
      , vdco{ (1.0 / std::pow( dy, 2) ) * KV }
    {
      for( int ix = 0; ix < NX; ix++ ) {
        for( int jy = 0; jy < NY; jy++ ) {
          P[ix][jy] = newDenseMat( NS, NS );
          Jbd[ix][jy] = newDenseMat( NS, NS );
          pivot[ix][jy] = newIndexArray( NS );
        }
      }
    }

    ~Data() { 
      for( int ix = 0; ix < NX; ix++ ) {
        for( int jy = 0; jy < NY; jy++ ) {
          destroyMat( P[ix][jy] );
          destroyMat( Jbd[ix][jy] );
          destroyArray( pivot[ix][jy] );
        }
      }
    }
  };

  static void SetInitialProfiles( std::vector<realtype>& state
                                  , realtype const dx
                                  , realtype const dy) {
    
    for( int iy = 0; iy < NY; iy++ ) {
      realtype const y = Y_MIN + iy * dy;
      realtype cy = std::pow( 0.1 * ( y - Y_MID ), 2 );
      cy = 1.0 - cy + 0.5 * std::pow( cy, 2 );
      for( int ix = 0; ix < NX; ix++ ) {
	realtype const x = X_MIN + ix * dx;
	realtype cx = std::pow( 0.1 * ( x - X_MID ), 2 );
	cx = 1.0 - cx + 0.5 * std::pow( cx, 2 );
	IJKth( state, 1, ix, iy) = SCALE_C_1 * cx * cy; 
	IJKth( state, 2, ix, iy) = SCALE_C_2 * cx * cy;
      }
    }
  }
  
  // right-hand-side routine
  static int RightHandSide( realtype time
                            , N_Vector state_vec
                            , N_Vector rhs_vec
                            , void *user_data ) {
    // get pointers to state (in)
    realtype const * const state  = N_VGetArrayPointer( state_vec );
    // get pointers to rhs (in/out)
    realtype* const rhs = N_VGetArrayPointer( rhs_vec );
    // get pointers to data (in/out)
    Data* const data = static_cast<Data* const>( user_data );

    // Set diurnal rate coefficients
    realtype q_3;
    {
      realtype const s = std::sin( data->om * time );
      if( s > 0.0 ) {
        q_3 = std::exp( -COEF_A_3 / s );
        data->q_4 = std::exp( -COEF_A_4 / s );
      } else {
        q_3 = 0.0;
        data->q_4 = 0.0;
      }
    }

    // Loop over all grid points

    for( int iy = 0; iy < NY; iy++ ) {

      // Set vertical diffusion coefficients at iy +- 1/2 

      realtype const y_dn = Y_MIN + ( iy - 0.5 ) * data->dy;
      realtype const y_up = y_dn + data->dy;
      realtype const cy_dn = data->vdco * std::exp( 0.2 * y_dn );
      realtype const cy_up = data->vdco * std::exp( 0.2 * y_up );
      
      int const iy_dn = iy + ( ( iy == 0    ) ? 1  : -1 );
      int const iy_up = iy + ( ( iy == NY-1 ) ? -1 : 1  );
      
      for( int ix = 0; ix < NX; ix++ ) {
  	// Extract c1 and c2, and set kinetic rate terms
  	realtype const c_1 = IJKth( state, 1, ix, iy ); 
  	realtype const c_2 = IJKth( state, 2, ix, iy );
  	realtype const qq_1 = COEF_Q_1 * c_1 * COEF_C_3;
  	realtype const qq_2 = COEF_Q_2 * c_1 * c_2;
  	realtype const qq_3 = q_3 * COEF_C_3;
  	realtype const qq_4 = data->q_4 * c_2;
  	realtype const rkin_1 = -qq_1 - qq_2 + 2.0 * qq_3 + qq_4;
  	realtype const rkin_2 = qq_1 - qq_2 - qq_4;

  	// Set vertical diffusion terms
  	realtype const c_1_dn = IJKth( state, 1, ix, iy_dn );
  	realtype const c_2_dn = IJKth( state, 2, ix, iy_dn );
  	realtype const c_1_up = IJKth( state, 1, ix, iy_up );
  	realtype const c_2_up = IJKth( state, 2, ix, iy_up );
  	realtype const vert_d_1 = cy_up * ( c_1_up - c_1 ) - cy_dn * ( c_1 - c_1_dn );
  	realtype const vert_d_2 = cy_up * ( c_2_up - c_2 ) - cy_dn * ( c_2 - c_2_dn );

  	// Set horizontal diffusion and advection terms
  	int const ix_left  = ix + ( ( ix == 0    ) ? 1  : -1 );
  	int const ix_right = ix + ( ( ix == NX-1 ) ? -1 : 1  );
  	realtype const c_1_lt = IJKth( state, 1, ix_left,  iy );
  	realtype const c_2_lt = IJKth( state, 2, ix_left,  iy );
  	realtype const c_1_rt = IJKth( state, 1, ix_right, iy );
  	realtype const c_2_rt = IJKth( state, 2, ix_right, iy );
  	realtype const hor_d_1 = data->hdco * ( c_1_rt - 2.0 * c_1 + c_1_lt );
  	realtype const hor_d_2 = data->hdco * ( c_2_rt - 2.0 * c_2 + c_2_lt );
  	realtype const hora_d_1 = data->haco * ( c_1_rt - c_1_lt );
  	realtype const hora_d_2 = data->haco * ( c_2_rt - c_2_lt );

  	// Load all terms into udot
  	IJKth( rhs, 1, ix, iy ) = vert_d_1 + hor_d_1 + hora_d_1 + rkin_1; 
  	IJKth( rhs, 2, ix, iy ) = vert_d_2 + hor_d_2 + hora_d_2 + rkin_2;
      }
    }

    return EXIT_SUCCESS;
  }

  // Jacobian-times-vector routine
  static int JacobianTimesVector( N_Vector vector_nvec
                                  , N_Vector jacobian_vector_nvec
                                  , realtype time
                                  , N_Vector state_nvec
                                  , N_Vector /*rhs_nvec*/
                                  , void *user_data
                                  , N_Vector /*tmp*/ ) {
    // get pointer to data (in/out)
    Data* const data = static_cast<Data* const>( user_data );

    realtype const * const state  = N_VGetArrayPointer( state_nvec );
    realtype const * const vector  = N_VGetArrayPointer( vector_nvec );
    realtype* const jacoabian_vector = N_VGetArrayPointer( jacobian_vector_nvec );

    // Set diurnal rate coefficients
    {
      double const s = std::sin( data->om * time );
      data->q_4 = ( s > 0.0 ) ? std::exp(-COEF_A_4 / s) : 0.0;
    }

    // Loop over all grid points

    for( int iy = 0; iy < NY; iy++ ) {

      // Set vertical diffusion coefficients at iy +- 1/2

      realtype const ydn = Y_MIN + (iy - 0.5) * data->dy;
      realtype const yup = ydn + data->dy;

      realtype const cydn = data->vdco * std::exp( 0.2 * ydn);
      realtype const cyup = data->vdco * std::exp( 0.2 * yup);

      int const iy_dn = iy + ( ( iy == 0    ) ? 1  : -1 );
      int const iy_up = iy + ( ( iy == NY-1 ) ? -1 : 1  );

      for( int ix = 0; ix < NX; ix++ ) {

        // Extract c1 and c2 at the current location and at neighbors

  	realtype const c1 = IJKth( state, 1, ix, iy ); 
  	realtype const c2 = IJKth( state, 2, ix, iy );

  	realtype const v1 = IJKth( vector, 1, ix, iy ); 
        realtype const v2 = IJKth( vector, 2, ix, iy );

  	realtype const c1dn = IJKth( state, 1, ix, iy_dn );
  	realtype const c2dn = IJKth( state, 2, ix, iy_dn );
  	realtype const c1up = IJKth( state, 1, ix, iy_up );
  	realtype const c2up = IJKth( state, 2, ix, iy_up );

  	realtype const v1dn = IJKth( vector, 1, ix, iy_dn );
  	realtype const v2dn = IJKth( vector, 2, ix, iy_dn );
  	realtype const v1up = IJKth( vector, 1, ix, iy_up );
  	realtype const v2up = IJKth( vector, 2, ix, iy_up );

  	int const ix_left  = ix + ( ( ix == 0    ) ? 1  : -1 );
  	int const ix_right = ix + ( ( ix == NX-1 ) ? -1 : 1  );

  	realtype const c1lt = IJKth( state, 1, ix_left,  iy ); 
  	realtype const c2lt = IJKth( state, 2, ix_left,  iy );
  	realtype const c1rt = IJKth( state, 1, ix_right, iy );
  	realtype const c2rt = IJKth( state, 2, ix_right, iy );

  	realtype const v1lt = IJKth( vector, 1, ix_left,  iy ); 
  	realtype const v2lt = IJKth( vector, 2, ix_left,  iy );
  	realtype const v1rt = IJKth( vector, 1, ix_right, iy );
  	realtype const v2rt = IJKth( vector, 2, ix_right, iy );

  	// Set kinetic rate terms

  	/* 
  	   rkin1 = -COEF_Q_1*COEF_C_3 * c1 - COEF_Q_2 * c1*c2 
           + data->q_4 * c2  + 2.0*COEF_C_3*q3;
  	   rkin2 =  COEF_Q_1*COEF_C_3 * c1 - COEF_Q_2 * c1*c2 - data->q_4 * c2; 
  	*/

        realtype Jv1 = 0.0;
  	realtype Jv2 = 0.0;
        
  	Jv1 += -( COEF_Q_1 * COEF_C_3 + COEF_Q_2 * c2 ) * v1
          + ( data->q_4 - COEF_Q_2 * c1 ) * v2;
  	Jv2 +=  ( COEF_Q_1 * COEF_C_3 - COEF_Q_2 * c2 ) * v1
          -  ( data->q_4 + COEF_Q_2 * c1 ) * v2;

  	// Set vertical diffusion terms

  	/* 
  	   vertd1 = -(cyup+cydn) * c1 + cyup * c1up + cydn * c1dn;
  	   vertd2 = -(cyup+cydn) * c2 + cyup * c2up + cydn * c2dn;
  	*/
  	Jv1 += -( cyup + cydn ) * v1  +  cyup * v1up  +  cydn * v1dn;
  	Jv2 += -( cyup + cydn ) * v2  +  cyup * v2up  +  cydn * v2dn;

  	// Set horizontal diffusion and advection terms

  	/* 
  	   hord1 = data->hdco*(c1rt - 2.0*c1 + c1lt);
  	   hord2 = data->hoco*(c2rt - 2.0*c2 + c2lt);
  	*/
  	Jv1 += data->hdco * ( v1rt - 2.0*v1 + v1lt );
  	Jv2 += data->hdco * ( v2rt - 2.0*v2 + v2lt );

  	/* 
  	   horad1 = data->haco*(c1rt - c1lt);
  	   horad2 = data->haco*(c2rt - c2lt);
  	*/
  	Jv1 += data->haco * ( v1rt - v1lt );
  	Jv2 += data->haco * ( v2rt - v2lt );

  	// Load two components of J*v

  	/* 
  	   IJKth(dstate, 1, ix, iy) = vertd1 + hord1 + horad1 + rkin1; 
  	   IJKth(dstate, 2, ix, iy) = vertd2 + hord2 + horad2 + rkin2;
  	*/
  	IJKth( jacoabian_vector, 1, ix, iy ) = Jv1;
  	IJKth( jacoabian_vector, 2, ix, iy ) = Jv2;

      }

    }

    return EXIT_SUCCESS;

  }

  // preconditioner setup routine
  static int PreconditionerSetup( realtype /*time*/
                                  , N_Vector state_nvec
                                  , N_Vector /*rhs_nvec*/
                                  , booleantype jac_ok
                                  , booleantype* jac_cur_ptr
                                  , realtype gamma
                                  , void* user_data) {
    // get pointer to data (in/out)
    Data* const data = static_cast<Data* const>( user_data );
    // get pointer to state (in)
    realtype const * const state = N_VGetArrayPointer( state_nvec);
    // extrat necessary variables from user data
    realtype **(*P)[NY] = data->P;
    realtype **(*Jbd)[NY] = data->Jbd;
    sunindextype *(*pivot)[NY] = data->pivot;
    
    if( jac_ok ) {

      // jac is OK: copy Jbd to P
      for( int iy = 0; iy < NY; iy++ ) {
        for( int ix = 0; ix < NX; ix++ ) {
          denseCopy( Jbd[ix][iy]
                     , P[ix][iy]
                     , NS
                     , NS);
        }
      }
      
      *jac_cur_ptr = SUNFALSE;

    } else {
      // jac needs to be update: generate Jbd from scratch and copy to P
      
      // Compute 2x2 diagonal Jacobian blocks (using q_4 values
      // computed on the last f call).  Load into P.
      for( int iy = 0; iy < NY; iy++ ) {
        realtype const y_dn = Y_MIN + (iy - 0.5 ) * data->dy;
        realtype const y_up = y_dn + data->dy;
        realtype const cy_dn = data->vdco * std::exp( 0.2 * y_dn );
        realtype const cy_up = data->vdco * std::exp( 0.2 * y_up );
        realtype const diag = -( cy_dn + cy_up + 2.0 * data->hdco );
        for( int ix = 0; ix < NX; ix++ ) {
          realtype const c1 = IJKth( state, 1, ix, iy );
          realtype const c2 = IJKth( state, 2, ix, iy );
          realtype** j = Jbd[ix][iy];
          realtype** a = P[ix][iy];
          j[0][0] = ( -COEF_Q_1 * COEF_C_3 - COEF_Q_2 * c2) + diag;
          j[1][0] = -COEF_Q_2 * c1 + data->q_4;
          j[0][1] = COEF_Q_1 * COEF_C_3 - COEF_Q_2 * c2;
          j[1][1] = ( -COEF_Q_2 * c1 - data->q_4 ) + diag;          
          denseCopy( j, a, NS, NS);
        }
      }

      *jac_cur_ptr = SUNTRUE;

    }

    // Scale by -gamma
    for( int iy = 0; iy < NY; iy++ ) {
      for( int ix = 0; ix < NX; ix++ ) {
        denseScale( -gamma, P[ix][iy], NS, NS);
      }
    }

    // Add identity matrix and do LU decompositions on blocks in place
    for( int ix = 0; ix < NX; ix++ ) {
      for( int iy = 0; iy < NY; iy++ ) {
        denseAddIdentity(P[ix][iy], NS);
        sunindextype retval = denseGETRF( P[ix][iy]
                                          , NS
                                          , NS
                                          , pivot[ix][iy] );
        if( retval != 0 ) return(1);
      }
    }
    
    return EXIT_SUCCESS;
  }
  
  static int PreconditionerSolve( realtype /*time*/
                                  , N_Vector /*state*/
                                  , N_Vector /*rhs*/
                                  , N_Vector r
                                  , N_Vector z
                                  , realtype /*gamma*/
                                  , realtype /*delta*/
                                  , int /*lr*/
                                  , void *user_data ) {    
    // get pointers to data (in/out)
    Data* const data = static_cast<Data* const>( user_data );
    
    // Extract the P and pivot arrays from user_data
    realtype **(*P)[NY] = data->P;
    sunindextype *(*pivot)[NY] = data->pivot;
    
    realtype * z_ptr = N_VGetArrayPointer( z );

    // copy r to z
    N_VScale( 1.0, r, z );
    
    // Solve the block-diagonal system Px = r using LU factors stored
    // in P and pivot data in pivot, and return the solution in z.
    
    for( int jx = 0; jx < NX; jx++ ) {
      for( int jy = 0; jy < NY; jy++ ) {
        realtype* v = &( IJKth( z_ptr, 1, jx, jy ) );
        denseGETRS( P[jx][jy]
                    , NS
                    , pivot[jx][jy]
                    , v );
      }
    }
    
    return EXIT_SUCCESS;
  }
};

TEST_CASE( "2-species diurnal kinetics advection-diffusion PDE "
           "system in 2 space dimensions can be solved", "[solver]" ) {
  
  SUNDIALS::SPGMR strategy;
  strategy.SetUserFucntions( nullptr
                             , CVDiurnalKry::JacobianTimesVector
                             , CVDiurnalKry::PreconditionerSetup
                             , CVDiurnalKry::PreconditionerSolve );
  strategy.SetLinearSolverOptions( SUNDIALS::Types::Preconditioner::left
                                   , SUNDIALS::Types::GramSchmidt::modified
                                   , 5
                                   , 0 );
      
  // set strategy
  SUNDIALS::CVODE cvode( &strategy );
      
  // set user functions
  cvode.SetUserFucntions( CVDiurnalKry::RightHandSide
                          , nullptr
                          , nullptr );
      
  // set linear multistep method
  cvode.SetLinearMultiStepMethod( SUNDIALS::Types::LinearMultisptepMethod::BDF );
      
  // set user data
  CVDiurnalKry::Data data;
  cvode.SetUserData( &data );
      
  // initial conditions
  std::vector<realtype> state( CVDiurnalKry::NEQ, 0.0 );
  CVDiurnalKry::SetInitialProfiles( state, data.dx, data.dy );
  cvode.SetInitialConditions( state );

  // set tolerances: optional
  cvode.SetTolerances( CVDiurnalKry::REL_TOL, CVDiurnalKry::ABS_TOL );

  // integration time
  cvode.SetIntegrationTime( 0.0, std::numeric_limits<double>::infinity() );

  // initialise
  REQUIRE( 0 == cvode.Initialise() );

  // integrate
  realtype t_out = 0.0;
  std::vector<realtype> test_state;
  for ( int i = 1; i <= CVDiurnalKry::NOUT; i++ ) {
    // set target time
    t_out += CVDiurnalKry::TWO_HR;
        
    // integrate
    REQUIRE( CV_SUCCESS == cvode.Integrate( t_out, state ) );
        
    // append the vector of states state
    test_state.insert( test_state.end()
                       , state.begin()
                       , state.end());

    // log if verbose
    if( TestConfig::verbose ) {
      int const left  = 0;
      int const right = CVDiurnalKry::NX - 1;
      int const bot   = 0;
      int const top   = CVDiurnalKry::NY - 1;
      int const mid_h = CVDiurnalKry::NX/2 - 1;
      int const mid_v = CVDiurnalKry::NY/2 - 1;
      std::cout << fmt::format( "t = {:.2e} \n"
                                "c1 (bot.left/middle/top rt.) = {:15.8e}  {:15.8e}  {:15.8e}\n"
                                "c2 (bot.left/middle/top rt.) = {:15.8e}  {:15.8e}  {:15.8e}\n\n"
                                , t_out
                                , IJKth( state, 1, left , bot   )
                                , IJKth( state, 1, mid_h, mid_v )
                                , IJKth( state, 1, right, top   )
                                , IJKth( state, 2, left , bot   )
                                , IJKth( state, 2, mid_h, mid_v )
                                , IJKth( state, 2, right, top   ) );

    }
  }
  
  if( TestConfig::verbose ) {
    std::cout << cvode.PrintSolverStatistics();
  }

  // read and store expected state from test data dir
  std::vector<realtype> expected_state;
  {
    std::ifstream bench_file;
    realtype val;
    bench_file.open( "../tests/data/cvDiurnal_kry.out", std::ifstream::in );
    if( bench_file.is_open() ) {
      while ( bench_file >> val )  {
        expected_state.push_back( val );
      }
      bench_file.close();
    }
  }

  // computed and expected state vectors must have the same size
  REQUIRE( test_state.size() == expected_state.size());

  // computed and expected state vectors must be equal within a tiny tolerance
  REQUIRE_THAT( test_state
                , Catch::Approx( expected_state ).epsilon( TestConfig::tolerance ) );

}
