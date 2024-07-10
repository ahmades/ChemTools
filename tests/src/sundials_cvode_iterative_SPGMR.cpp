#include <algorithm>
#include <fstream>
#include <limits>
#include <vector>

#include <cmath>

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/printf.h>

#include <catch2/catch.hpp>

#include "sundials/cvode/interface.hpp"
#include "test_config.hpp"

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

struct Dimension
{
  realtype min;
  realtype max;
  realtype mid;
  int npts;
  realtype delta;

  Dimension(realtype min_, realtype max_, int npts_)
      : min{min_}, max{max_}, mid{(min + max) / 2.0}, npts{npts_}, delta{(max_ - min_) / (npts_ - 1)}
  {
  }
};

struct Domain
{
  Dimension dim_x;
  Dimension dim_y;

  Domain(Dimension dim_x_, Dimension dim_y_)
      : dim_x(dim_x_), dim_y(dim_y_)
  {
  }
};

struct PhysicalParameters
{
  realtype Kh;            // horizontal diffusivity Kh
  realtype coef_Kv;       // coefficient in the vertical diffusivity, Kv
  realtype velocity;      // advection velocity V
  realtype coef_q_1;      // coefficient q1
  realtype coef_q_2;      // coefficient q2
  realtype coef_c_3;      // coefficient c3
  realtype coef_a_3;      // coefficient in expression for q3(t)
  realtype coef_a_4;      // coefficient in expression for q4(t)
  realtype coef_init_c_1; // coefficients 1 in initial profiles
  realtype coef_init_c_2; // coefficients 2 in initial profiles

  PhysicalParameters(realtype Kh_, realtype coef_Kv_, realtype velocity_, realtype coef_q_1_, realtype coef_q_2_, realtype coef_c_3_, realtype coef_a_3_, realtype coef_a_4_, realtype coef_init_c_1_, realtype coef_init_c_2_)
      : Kh{Kh_}, coef_Kv{coef_Kv_}, velocity{velocity_}, coef_q_1{coef_q_1_}, coef_q_2{coef_q_2_}, coef_c_3{coef_c_3_}, coef_a_3{coef_a_3_}, coef_a_4{coef_a_4_}, coef_init_c_1{coef_init_c_1_}, coef_init_c_2{coef_init_c_2_}
  {
  }
};

struct Time
{
  realtype start;
  realtype stop;
  Time(realtype start_, realtype stop_)
      : start(start_), stop(stop_)
  {
  }
};

class CVDiurnalKry : public SUNDIALS::CVODE::Client
{

public:
  CVDiurnalKry(Domain const& domain, PhysicalParameters const& phys_params, Time const& time, SUNDIALS::CVODE::Types::IntegrationTolerance const& tolerance)
      : m_domain(domain), m_phys_params(phys_params), m_time(time), m_tolerance(tolerance), m_neq{m_nspecies * m_domain.dim_x.npts * m_domain.dim_y.npts}, m_coef_q_4{0.0}, m_frequency{PI / (12.0 * 3600.0)}, m_hdco{phys_params.Kh / std::pow(m_domain.dim_x.delta, 2)}, m_vdco{(1.0 / std::pow(m_domain.dim_y.delta, 2)) * m_phys_params.coef_Kv}, m_haco{phys_params.velocity / (2.0 * m_domain.dim_x.delta)}, m_prec_mat(m_domain.dim_x.npts, std::vector<realtype**>(m_domain.dim_y.npts)), m_jac_mat(m_domain.dim_x.npts, std::vector<realtype**>(m_domain.dim_y.npts)), m_pivot(m_domain.dim_x.npts, std::vector<sunindextype*>(m_domain.dim_y.npts))
  {
    for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
      {
        for (int jy = 0; jy < m_domain.dim_y.npts; jy++)
          {
            m_prec_mat[ix][jy] = newDenseMat(m_nspecies, m_nspecies);
            m_jac_mat[ix][jy] = newDenseMat(m_nspecies, m_nspecies);
            m_pivot[ix][jy] = newIndexArray(m_nspecies);
          }
      }

    // enable optional user functions
    Client::opt_udf_set.jacobian_times_vector = true;
  }

  ~CVDiurnalKry()
  {
    for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
      {
        for (int jy = 0; jy < m_domain.dim_y.npts; jy++)
          {
            destroyMat(m_prec_mat[ix][jy]);
            destroyMat(m_jac_mat[ix][jy]);
            destroyArray(m_pivot[ix][jy]);
          }
      }
  }

  void Initialise()
  {
    // set initial conditions
    {
      std::vector<realtype> state(m_neq, 0.0);
      for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
        {
          realtype const y = m_domain.dim_y.min + iy * m_domain.dim_y.delta;
          realtype cy = std::pow(0.1 * (y - m_domain.dim_y.mid), 2);
          cy = 1.0 - cy + 0.5 * std::pow(cy, 2);
          for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
            {
              realtype const x = m_domain.dim_x.min + ix * m_domain.dim_x.delta;
              realtype cx = std::pow(0.1 * (x - m_domain.dim_x.mid), 2);
              cx = 1.0 - cx + 0.5 * std::pow(cx, 2);
              AccessData(state, 1, ix, iy) = m_phys_params.coef_init_c_1 * cx * cy;
              AccessData(state, 2, ix, iy) = m_phys_params.coef_init_c_2 * cx * cy;
            }
        }
      SetState(state);
    }

    // set tolerances: optional
    SetTolerances(m_tolerance.relative, m_tolerance.absolute);

    // integration time
    SetIntegrationTime(m_time.start, m_time.stop);
  }

  // very dangerous
  realtype& AccessData(realtype* data, int i, int j, int k)
  {
    assert(i >= 1 && i <= m_nspecies);
    assert(j >= 0 && j < m_domain.dim_x.npts);
    assert(k >= 0 && k < m_domain.dim_y.npts);
    return data[i - 1
                + j * m_nspecies
                + k * m_nspecies * m_domain.dim_x.npts];
  }

  realtype& AccessData(std::vector<realtype>& data, int i, int j, int k)
  {
    return AccessData(data.data(), i, j, k);
  }

private:
  static int constexpr m_nspecies = 2;            // number of species
  static realtype constexpr PI = 3.1415926535898; // not using M_PI from cmath to match sundial's results

  Domain m_domain;                                          // input: domain parameters
  PhysicalParameters m_phys_params;                         // input: physical parameters
  Time m_time;                                              // input: time info
  SUNDIALS::CVODE::Types::IntegrationTolerance m_tolerance; // input: tolerances

  int m_neq;                                       // number of equations
  realtype m_coef_q_4;                             // coefficient q4
  realtype m_frequency;                            // angular frequency in thr rate constant
  realtype m_hdco;                                 // horizontal diffusion discretisation coef
  realtype m_vdco;                                 // vertical diffusion discretisation coef
  realtype m_haco;                                 // horizontal advection discretisation coef
  std::vector<std::vector<realtype**>> m_prec_mat; // preconditioning matrices
  std::vector<std::vector<realtype**>> m_jac_mat;  // jacobain matrices
  std::vector<std::vector<sunindextype*>> m_pivot; // pivot vectors

  // right-hand-side routine
  int RightHandSide(realtype const time, realtype* const state, realtype* rhs) override
  {
    // Set diurnal rate coefficients
    realtype coef_q_3;
    {
      realtype const s = std::sin(m_frequency * time);
      if (s > 0.0)
        {
          coef_q_3 = std::exp(-m_phys_params.coef_a_3 / s);
          m_coef_q_4 = std::exp(-m_phys_params.coef_a_4 / s);
        }
      else
        {
          coef_q_3 = 0.0;
          m_coef_q_4 = 0.0;
        }
    }

    // Loop over all grid points

    for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
      {

        // Set vertical diffusion coefficients at iy +- 1/2

        realtype const y_dn = m_domain.dim_y.min + (iy - 0.5) * m_domain.dim_y.delta;
        realtype const y_up = y_dn + m_domain.dim_y.delta;
        realtype const cy_dn = m_vdco * std::exp(0.2 * y_dn);
        realtype const cy_up = m_vdco * std::exp(0.2 * y_up);

        int const iy_dn = iy + ((iy == 0) ? 1 : -1);
        int const iy_up = iy + ((iy == m_domain.dim_y.npts - 1) ? -1 : 1);

        for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
          {
            // Extract c1 and c2, and set kinetic rate terms
            realtype const c_1 = AccessData(state, 1, ix, iy);
            realtype const c_2 = AccessData(state, 2, ix, iy);

            realtype const qq_1 = m_phys_params.coef_q_1 * c_1 * m_phys_params.coef_c_3;
            realtype const qq_2 = m_phys_params.coef_q_2 * c_1 * c_2;
            realtype const qq_3 = coef_q_3 * m_phys_params.coef_c_3;
            realtype const qq_4 = m_coef_q_4 * c_2;

            realtype const rkin_1 = -qq_1 - qq_2 + 2.0 * qq_3 + qq_4;
            realtype const rkin_2 = qq_1 - qq_2 - qq_4;

            // Set vertical diffusion terms
            realtype const c_1_dn = AccessData(state, 1, ix, iy_dn);
            realtype const c_2_dn = AccessData(state, 2, ix, iy_dn);

            realtype const c_1_up = AccessData(state, 1, ix, iy_up);
            realtype const c_2_up = AccessData(state, 2, ix, iy_up);

            realtype const vert_d_1 = cy_up * (c_1_up - c_1) - cy_dn * (c_1 - c_1_dn);
            realtype const vert_d_2 = cy_up * (c_2_up - c_2) - cy_dn * (c_2 - c_2_dn);

            // Set horizontal diffusion and advection terms
            int const ix_left = ix + ((ix == 0) ? 1 : -1);
            int const ix_right = ix + ((ix == m_domain.dim_x.npts - 1) ? -1 : 1);

            realtype const c_1_lt = AccessData(state, 1, ix_left, iy);
            realtype const c_2_lt = AccessData(state, 2, ix_left, iy);
            realtype const c_1_rt = AccessData(state, 1, ix_right, iy);
            realtype const c_2_rt = AccessData(state, 2, ix_right, iy);

            realtype const hor_d_1 = m_hdco * (c_1_rt - 2.0 * c_1 + c_1_lt);
            realtype const hor_d_2 = m_hdco * (c_2_rt - 2.0 * c_2 + c_2_lt);
            realtype const hora_d_1 = m_haco * (c_1_rt - c_1_lt);
            realtype const hora_d_2 = m_haco * (c_2_rt - c_2_lt);

            // Load all terms into udot
            AccessData(rhs, 1, ix, iy) = vert_d_1 + hor_d_1 + hora_d_1 + rkin_1;
            AccessData(rhs, 2, ix, iy) = vert_d_2 + hor_d_2 + hora_d_2 + rkin_2;
          }
      }

    return EXIT_SUCCESS;
  }

  // Jacobian-times-vector routine
  int JacobianTimesVector(realtype* const vector, realtype* jacobian_vector, realtype const time, realtype* const state, realtype* const /*rhs_nvec*/) override
  {

    // Set diurnal rate coefficients
    {
      double const s = std::sin(m_frequency * time);
      m_coef_q_4 = (s > 0.0) ? std::exp(-m_phys_params.coef_a_4 / s) : 0.0;
    }

    // Loop over all grid points

    for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
      {

        // Set vertical diffusion coefficients at iy +- 1/2

        realtype const ydn = m_domain.dim_y.min + (iy - 0.5) * m_domain.dim_y.delta;
        realtype const yup = ydn + m_domain.dim_y.delta;

        realtype const cydn = m_vdco * std::exp(0.2 * ydn);
        realtype const cyup = m_vdco * std::exp(0.2 * yup);

        int const iy_dn = iy + ((iy == 0) ? 1 : -1);
        int const iy_up = iy + ((iy == m_domain.dim_y.npts - 1) ? -1 : 1);

        for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
          {

            // Extract c1 and c2 at the current location and at neighbors

            realtype const c1 = AccessData(state, 1, ix, iy);
            realtype const c2 = AccessData(state, 2, ix, iy);

            realtype const v1 = AccessData(vector, 1, ix, iy);
            realtype const v2 = AccessData(vector, 2, ix, iy);

            realtype const v1dn = AccessData(vector, 1, ix, iy_dn);
            realtype const v2dn = AccessData(vector, 2, ix, iy_dn);

            realtype const v1up = AccessData(vector, 1, ix, iy_up);
            realtype const v2up = AccessData(vector, 2, ix, iy_up);

            int const ix_left = ix + ((ix == 0) ? 1 : -1);
            int const ix_right = ix + ((ix == m_domain.dim_x.npts - 1) ? -1 : 1);

            realtype const v1lt = AccessData(vector, 1, ix_left, iy);
            realtype const v2lt = AccessData(vector, 2, ix_left, iy);

            realtype const v1rt = AccessData(vector, 1, ix_right, iy);
            realtype const v2rt = AccessData(vector, 2, ix_right, iy);

            // Set kinetic rate terms

            /*
               rkin1 = -m_phys_params.coef_q_1*m_phys_params.coef_c_3 * c1
               - m_phys_params.coef_q_2 * c1*c2
               + m_coef_q_4 * c2  + 2.0*m_phys_params.coef_c_3*q3;
               rkin2 =  m_phys_params.coef_q_1*m_phys_params.coef_c_3 * c1
               - m_phys_params.coef_q_2 * c1*c2 - m_coef_q_4 * c2;
            */

            realtype Jv1 = 0.0;
            realtype Jv2 = 0.0;

            Jv1 += -(m_phys_params.coef_q_1 * m_phys_params.coef_c_3
                     + m_phys_params.coef_q_2 * c2)
                       * v1
                   + (m_coef_q_4 - m_phys_params.coef_q_2 * c1) * v2;
            Jv2 += (m_phys_params.coef_q_1 * m_phys_params.coef_c_3
                    - m_phys_params.coef_q_2 * c2)
                       * v1
                   - (m_coef_q_4 + m_phys_params.coef_q_2 * c1) * v2;

            // Set vertical diffusion terms

            /*
               vertd1 = -(cyup+cydn) * c1 + cyup * c1up + cydn * c1dn;
               vertd2 = -(cyup+cydn) * c2 + cyup * c2up + cydn * c2dn;
            */
            Jv1 += -(cyup + cydn) * v1 + cyup * v1up + cydn * v1dn;
            Jv2 += -(cyup + cydn) * v2 + cyup * v2up + cydn * v2dn;

            // Set horizontal diffusion and advection terms

            /*
               hord1 = m_hdco*(c1rt - 2.0*c1 + c1lt);
               hord2 = m_hoco*(c2rt - 2.0*c2 + c2lt);
            */
            Jv1 += m_hdco * (v1rt - 2.0 * v1 + v1lt);
            Jv2 += m_hdco * (v2rt - 2.0 * v2 + v2lt);

            /*
               horad1 = m_haco*(c1rt - c1lt);
               horad2 = m_haco*(c2rt - c2lt);
            */
            Jv1 += m_haco * (v1rt - v1lt);
            Jv2 += m_haco * (v2rt - v2lt);

            // Load two components of J*v

            /*
               AccessData(dstate, 1, ix, iy) = vertd1 + hord1 + horad1 + rkin1;
               AccessData(dstate, 2, ix, iy) = vertd2 + hord2 + horad2 + rkin2;
            */
            AccessData(jacobian_vector, 1, ix, iy) = Jv1;
            AccessData(jacobian_vector, 2, ix, iy) = Jv2;
          }
      }

    return EXIT_SUCCESS;
  }

  // preconditioner setup routine
  int PreconditionerSetup(realtype const /*time*/
                          ,
                          realtype* const state,
                          realtype* const /*rhs*/
                          ,
                          booleantype const jac_ok,
                          booleantype* jac_cur_ptr,
                          realtype const gamma)
  {
    if (jac_ok)
      {

        // jac is OK: copy Jbd to P
        for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
          {
            for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
              {
                denseCopy(m_jac_mat[ix][iy], m_prec_mat[ix][iy], m_nspecies, m_nspecies);
              }
          }

        *jac_cur_ptr = SUNFALSE;
      }
    else
      {
        // jac needs to be update: generate Jbd from scratch and copy to P

        // Compute 2x2 diagonal Jacobian blocks (using q_4 values
        // computed on the last f call).  Load into P.
        for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
          {
            realtype const y_dn = m_domain.dim_y.min + (iy - 0.5) * m_domain.dim_y.delta;
            realtype const y_up = y_dn + m_domain.dim_y.delta;
            realtype const cy_dn = m_vdco * std::exp(0.2 * y_dn);
            realtype const cy_up = m_vdco * std::exp(0.2 * y_up);
            realtype const diag = -(cy_dn + cy_up + 2.0 * m_hdco);
            for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
              {
                realtype const c1 = AccessData(state, 1, ix, iy);
                realtype const c2 = AccessData(state, 2, ix, iy);
                realtype** j = m_jac_mat[ix][iy];
                realtype** a = m_prec_mat[ix][iy];
                j[0][0] = (-m_phys_params.coef_q_1 * m_phys_params.coef_c_3
                           - m_phys_params.coef_q_2 * c2)
                          + diag;
                j[1][0] = -m_phys_params.coef_q_2 * c1 + m_coef_q_4;
                j[0][1] = m_phys_params.coef_q_1 * m_phys_params.coef_c_3 - m_phys_params.coef_q_2 * c2;
                j[1][1] = (-m_phys_params.coef_q_2 * c1 - m_coef_q_4) + diag;
                denseCopy(j, a, m_nspecies, m_nspecies);
              }
          }

        *jac_cur_ptr = SUNTRUE;
      }

    // Scale by -gamma
    for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
      {
        for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
          {
            denseScale(-gamma, m_prec_mat[ix][iy], m_nspecies, m_nspecies);
          }
      }

    // Add identity matrix and do LU decompositions on blocks in place
    for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
      {
        for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
          {
            denseAddIdentity(m_prec_mat[ix][iy], m_nspecies);
            sunindextype retval = denseGETRF(m_prec_mat[ix][iy], m_nspecies, m_nspecies, m_pivot[ix][iy]);
            if (retval != 0)
              return (1);
          }
      }

    return EXIT_SUCCESS;
  }

  int PreconditionerSolve(realtype const /*time*/
                          ,
                          realtype* const /*state*/
                          ,
                          realtype* const /*rhs*/
                          ,
                          realtype* const r,
                          realtype* z,
                          realtype const /*gamma*/
                          ,
                          realtype const /*delta*/
                          ,
                          int const /*lr*/)
  {
    // copy r to z
    // N_VScale( 1.0, r, z );
    std::copy(r, r + m_neq, z);

    // Solve the block-diagonal system Px = r using LU factors stored
    // in P and pivot data in pivot, and return the solution in z.

    for (int ix = 0; ix < m_domain.dim_x.npts; ix++)
      {
        for (int iy = 0; iy < m_domain.dim_y.npts; iy++)
          {
            realtype* v = &(AccessData(z, 1, ix, iy));
            denseGETRS(m_prec_mat[ix][iy], m_nspecies, m_pivot[ix][iy], v);
          }
      }

    return EXIT_SUCCESS;
  }
};

TEST_CASE("2-species diurnal kinetics advection-diffusion PDE "
          "system in 2 space dimensions can be solved",
          "[iterative solver]")
{

  // create instance of the problem to solve
  Domain const domain({0.0, 20.0, 10},   // dim_x: min, max, npts
                      {30.0, 50.0, 10}); // dim_y: min, max, npts

  PhysicalParameters const phys_params(4.0e-6 // horizontal diffusivity Kh
                                       ,
                                       1.0e-8 // coefficient in the vertical diffusivity Kv
                                       ,
                                       0.001 // advection velocity V
                                       ,
                                       1.63e-16 // coefficient q1
                                       ,
                                       4.66e-16 // coefficient q2
                                       ,
                                       3.7e+16 // coefficient c3
                                       ,
                                       22.62 // coefficient in expression for q3(t)
                                       ,
                                       7.601 // coefficient in expression for q4(t)
                                       ,
                                       1.0e+6 // coefficients 1 in initial profiles
                                       ,
                                       1.0e+12 // coefficients 1 in initial profiles
  );

  Time const time(0.0, std::numeric_limits<double>::infinity());

  realtype const rel_tol = 1.0e-5;          // scalar relative tolerance
  realtype const abs_tol = 100.0 * rel_tol; // scalar absolute tolerance
  SUNDIALS::CVODE::Types::IntegrationTolerance const tolerance(rel_tol, abs_tol);

  CVDiurnalKry client(domain, phys_params, time, tolerance);

  // set initial conditions (and constraints if applicable), integration time, and tolerances
  // can be even moved to body of CVDiurnalKry constructor
  client.Initialise();

  // select strategy
  SUNDIALS::CVODE::SPGMR strategy;

  strategy.SetLinearSolverOptions(SUNDIALS::CVODE::Types::Preconditioner::left, SUNDIALS::CVODE::Types::GramSchmidt::modified, 5, 0);

  // set strategy
  SUNDIALS::CVODE::Solver cvode(&strategy);

  // set linear multistep method
  cvode.SetLinearMultiStepMethod(SUNDIALS::CVODE::Types::LinearMultisptepMethod::BDF);

  // set client data
  cvode.SetClientData(&client);

  // set time step control options
  cvode.SetInitStepSize(0.0); // 0 sets default
  cvode.SetMinStepSize(0.0);  // 0 sets default
  cvode.SetMaxStepSize(0.0);  // 0 sets default
  cvode.SetMaxNumSteps(600);  // 0 sets default, 500. 600 > 500 => should not influence the test

  // set solver control options
  cvode.SetMaxOrder(5);         // same value as target, otherwise test will fail
  cvode.SetMaxWarnMessages(15); // target uses default, 10, does not influence the test
  cvode.SetStabilityLimitDetection(true);
  cvode.SetMaxErrorTestFailures(10);          // target uses default, 7, but test should pass
  cvode.SetMaxNonlinearIterations(4);         // target uses default, 3, but test should pass
  cvode.SetMaxConvergenceFailures(15);        // target uses default, 10, but test should pass
  cvode.SetNonlinConvergenceCoefficient(0.1); // same value as target, otherwise test will fail

  // initialise
  REQUIRE(0 == cvode.Initialise());

  // integrate
  {
    realtype time = 0.0;
    realtype const time_step = 2.0 * 3600.0;  // 2 hrs
    realtype const time_stop = 24.0 * 3600.0; // 1 day
    std::vector<realtype> test_state;         // vector for storage of temporal solutions

    while (time < time_stop)
      {
        // set target time
        time += time_step;

        // integrate
        REQUIRE(CV_SUCCESS == cvode.Integrate(time));

        // append the vector of states
        test_state.insert(test_state.end(), client.State(), client.State() + client.NStates());

        // log if verbose
        if (TestConfig::verbose)
          {
            int const left = 0;
            int const right = domain.dim_x.npts - 1;
            int const bot = 0;
            int const top = domain.dim_y.npts - 1;
            int const mid_h = domain.dim_x.npts / 2 - 1;
            int const mid_v = domain.dim_y.npts / 2 - 1;
            fmt::print("t = {:.2e} \n"
                       "c1 (bot.left/middle/top rt.) = {:15.8e}  {:15.8e}  {:15.8e}\n"
                       "c2 (bot.left/middle/top rt.) = {:15.8e}  {:15.8e}  {:15.8e}\n\n",
                       time,
                       client.AccessData(client.State(), 1, left, bot),
                       client.AccessData(client.State(), 1, mid_h, mid_v),
                       client.AccessData(client.State(), 1, right, top),
                       client.AccessData(client.State(), 2, left, bot),
                       client.AccessData(client.State(), 2, mid_h, mid_v),
                       client.AccessData(client.State(), 2, right, top));
          }
      }

    if (TestConfig::verbose)
      {
        std::string const title = fmt::format(fmt::emphasis::underline
                                                  | fmt::emphasis::bold,
                                              "Final statistics:");
        fmt::print("{}\n{}\n", title, cvode.PrintSolverStatistics());
      }

    // read and store expected state from test data dir
    std::vector<realtype> expected_state;
    {
      std::ifstream file;
      realtype val;
      file.open("../tests/data/sundials_cvode/iterative_spgmr_cvDiurnal_kry.dat", std::ifstream::in);
      if (file.is_open())
        {
          while (file >> val)
            {
              expected_state.push_back(val);
            }
          file.close();
        }
    }

    // computed and expected state vectors must have the same size
    REQUIRE(test_state.size() == expected_state.size());

    // computed and expected state vectors must be equal within a tiny tolerance
    REQUIRE_THAT(test_state, Catch::Approx(expected_state).epsilon(TestConfig::tolerance));
  }
}
