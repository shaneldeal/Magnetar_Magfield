//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file parker_my_IC.cpp
//  \brief Problem generator for Parker wind.
//


#define COMP_DT 1    // Compute and save to uov
#define COMP_DIV_B 1
#define COMP_SRC    1
// C/C++ headers
#include <algorithm>  // min, max
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>   // std::numeric_limits<float>::epsilon()
#include <sstream>
#include <string>
// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"                  // Globals
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../scalars/scalars.hpp"

// Configuration checking

// Static Variables
// critical radius, inner radius, central gravitational parameter,
// isothermal sound speed, initial phi-velocity, initial field strength,
// tilt angle, frame angular velocity
static Real r_0, inv_r2, rho_0, mu, Ye, Na, p_0, Tc, T_0, dpdd_0, B_0, alpha,a_iso, GuessTemp_Bc,GuessYe_Bc, r0sq, fixed_tempBC;
static Real Coeff_nu_0, Coeff_nubar_0,Coeff_nu_f,Coeff_nubar_f,t_L_0, t_L_1, dCoeff_nu,L_nu,L_nubar,L_mu,L_nu_f,L_nubar_f,L_mu_f,eps_nu,eps_mu,eps_nubar,eps_nu_f,eps_mu_f,eps_nubar_f,dCoeff_nubar,Omega_0, energy_floor, Mpns, Omega_dot, r_env1,r_env2,env_mass,env_vel;
//static const Real float_eps = std::numeric_limits<Real>::epsilon();
static int rows; //number of rows in the IC file
static int IDT1, IDT2, IDT3, IDT4, IDT5;
static bool use_IC_file, set_T_at_r0,rot_flag,cs,env;
static int rot_option, source_choice, jdot_start, vr_choice, vth_choice, ye_index, temp_index;
static int ISRC1, ISRC2, ISRC3, ISRC4, ISRC5, ISRC6, ISRC7, ISRC8, ISRC9;
static int n_bc_ranks;
static MPI_Comm bc_comm;
  
void UserSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                    const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                    AthenaArray<Real> &cons_scalar);

Real grid_space_X2(Real x, RegionSize rs);

// Vector Potential
static Real A3(const Real x1, const Real x2, const Real x3) {
  Real a3 = 0.5 * B_0 * r_0 * std::pow(r_0/x1,2) *
    (std::sin(x2)*std::cos(alpha) - std::cos(x2)*std::cos(x3)*std::sin(alpha));
  return a3;
}

static Real A2(const Real x1, const Real x2, const Real x3) {
  Real a2 = -0.5 * B_0 * r_0 * std::pow(r_0/x1,2) * std::sin(x3) * std::sin(alpha);
  return a2;
}

static Real A1(const Real x1, const Real x2, const Real x3) {
  Real a1=0.0;
  return a1;
}


Real fermi_approx(int n, Real eta){
  if (n==2) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,3.0) / 3.0 + 3.2898681 * eta;
    Real ff = 2.0 * (a - std::pow(a,2.0) / 8.0 + std::pow(a,3.0) / 27.0 - std::pow(a,4.0) / 64.0 + std::pow(a,5.0) / 125.0 - std::pow(a,6.0) / 216.0);
    if (eta<0) {
        return ff;
    } else {
        return s+ff;
    }
  }
 else if (n==3) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,4.0) / 4.0 + 4.9348022 * std::pow(eta,2.0) + 11.351273;
    Real ff = 6.0 * (a - std::pow(a,2.0) / 16.0 + std::pow(a,3.0) / 81.0 - std::pow(a,4.0) / 256.0);
    if (eta<0) {
        return ff;
    } else {
        return s-ff;
    }
  }
  if (n==4) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,5.0) / 5.0 + 6.5797363 * std::pow(eta,3.0) + 45.457576 * eta;
    Real ff = 24.0 * (a - std::pow(a,2.0) / 32.0 + std::pow(a,3.0) / 243.0);
    if (eta<0.0) {
        return ff;
    } else {
        return s+ff;
    }
  }
  else if (n==5) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,6.0) / 6.0 + 8.2246703 * std::pow(eta,4.0) + 113.64394 * std::pow(eta,2.0) + 236.53226;
    Real ff = 120.0 * (a - std::pow(a,2.0) / 64.0 + std::pow(a,3.0) / 729.0);
    if(eta<0.0) {
        return ff;
    } else {
        return s-ff;
    }
  }
  else {
   std::stringstream msg;
    msg << "### FATAL ERROR in Mesh" << std::endl
        << "Argument different from 2, 3, 4 or 5 passed to Fermi function"<< std::endl;
    ATHENA_ERROR(msg);
    return 0.0;
  }
}


void lum_func(Real time, AthenaArray<Real> &out) {

  Real Coeff_nu;
  Real Coeff_nubar;
//Luminosity and epsilon
  Real lnu_y;
  Real epsnu_y;
  Real lnubar_y;
  Real epsnubar_y;
  Real lmu_y;
  Real epsmu_y;
  if (time<t_L_0) {
      lnu_y=L_nu;
      lmu_y=L_mu;
      lnubar_y=L_nubar;
      epsnu_y=eps_nu;
      epsmu_y=eps_mu;
      epsnubar_y=eps_nubar;
  }
  if (time>=t_L_0 && time<=t_L_1) {
      Real ti=time+0.7;
      lnu_y = L_nu*std::pow(ti/1.2,-0.982);
      lnubar_y = L_nubar*std::pow(ti/1.2,-0.969);
      lmu_y = L_mu*std::pow(ti/1.2,-0.8);  
    
      if(ti<=4.7655) {

        epsnubar_y=eps_nubar*std::pow(ti/1.2,-0.093);
        epsnu_y=eps_nu*std::pow(ti/1.2,-0.032);
        epsmu_y=eps_mu*std::pow(ti/1.2,-0.091);


      }
      else {
        Real f2=4.7655/1.2;
        epsnubar_y=eps_nubar*std::pow(f2,-0.093)*std::pow(lnubar_y/(L_nubar*std::pow(f2,-0.969)),0.25);
        epsnu_y=eps_nu*std::pow(f2,-0.032)*std::pow(lnu_y/(L_nu*std::pow(f2,-0.982)),0.25);
        epsmu_y=eps_mu*std::pow(f2,-0.091)*std::pow(lmu_y/(L_mu*std::pow(f2,-0.8)),0.25);
      }
  }
  if (time>t_L_1) {
      lnu_y=L_nu_f;
      lmu_y=L_mu_f;
      lnubar_y=L_nubar_f;
      epsnu_y=eps_nu_f;
      epsmu_y=eps_mu_f;
      epsnubar_y=eps_nubar_f;
  }


  Real ferm6 = 714.668;
  Real ferm5 = 118.266;
  Real ferm4 = 23.3309;
  Real ferm3 = 5.6822;
  Real ferm2 = 1.80309;
  
  Real ktemp_nubar = epsnubar_y*ferm2/ferm3;
  Real ktemp_nu = epsnu_y*ferm2/ferm3;


  Real eavg1_nubar = epsnubar_y;
  Real eavg2_nubar = ktemp_nubar*ktemp_nubar*ferm4/ferm2;
  

  Real eavg1_nu = epsnu_y;
  Real eavg2_nu = ktemp_nu*ktemp_nu*ferm4/ferm2;

  Coeff_nu = lnu_y*SQR(ktemp_nu)*ferm5/ferm3;
  Coeff_nubar = lnubar_y*SQR(ktemp_nubar)*ferm5/ferm3;

  out(0)=lnu_y;
  out(1)=lnubar_y;
  out(2)=lmu_y;
  out(3)=std::pow(ktemp_nu*ktemp_nu*ferm5/ferm3,0.5);
  out(4)=std::pow(ktemp_nubar*ktemp_nubar*ferm5/ferm3,0.5);
  out(5)=epsmu_y;
  out(6)=eavg1_nu;
  out(7)=eavg2_nu;
  out(8)=eavg1_nubar;
  out(9)=eavg2_nubar;
  out(10)=Coeff_nu;
  out(11)=Coeff_nubar;
}


// Heating/cooling function
Real qdotQW(Real temp, Real rho, Real x, Real ye, Real time, Real etaqw, AthenaArray<Real> &le, int return_index) {
  //return_index 0 for main and 1 for BCs
  // temp is in MeV
  // returns heating - cooling in units of MeV s^{-1} g^{-1} units


//Here epsnu_y and epsnubar_y are square roots of the energy term in charged current heating rate
  Real lnu_y = le(0);
  Real lnubar_y = le(1);
  Real lmu_y = le(2);
  Real epsnu_y = le(3);
  Real epsnubar_y = le(4);
  Real eavg1_mu = le(5);
  Real eavg1_nu = le(6) ;
  Real eavg2_nu = le(7);
  Real eavg1_nubar = le(8);
  Real eavg2_nubar = le(9);
  Real Coeff_nu = le(10);
  Real Coeff_nubar = le(11);

  Real ferm6 = 714.668;
  Real ferm5 = 118.266;
  Real ferm4 = 23.3309;
  Real ferm3 = 5.6822;
  Real ferm2 = 1.80309;
    

  Real h1=1e12*9.65*Na*((1.0-ye)*Coeff_nu + ye*Coeff_nubar)*(1.0-x)*inv_r2;

//neutrino electron positron scattering
  Real ferm2_eta = fermi_approx(2,etaqw);
  Real ferm3_eta = fermi_approx(3,etaqw);
  Real ferm4_eta = fermi_approx(4,etaqw);
  Real ferm2_negeta = fermi_approx(2,-1.0*etaqw);
  Real ferm3_negeta = fermi_approx(3,-1.0*etaqw);
  Real ferm4_negeta = fermi_approx(4,-1.0*etaqw);

  Real f1=1e8*(1.511/rho)*Na*std::pow(temp,4.0);
  Real h2_nubar = f1*lnubar_y*(0.926*ferm2_eta*(eavg1_nubar*ferm2/ferm3 - 4.0*temp*ferm3_eta/ferm4_eta) + 2.20858*ferm2_negeta*(eavg1_nubar*ferm2/ferm3 - 4.0*temp*ferm3_negeta/ferm4_negeta));

  Real h2_nu = f1*lnu_y*(2.20858*ferm2_eta*(eavg1_nu*ferm2/ferm3 - 4.0*temp*ferm3_eta/ferm4_eta) + 0.926*ferm2_negeta*(eavg1_nu*ferm2/ferm3 - 4.0*temp*ferm3_negeta/ferm4_negeta));
//factor of 2 for mu and tau
  Real h2_mu = 2.0*f1*lmu_y*(ferm2_eta*(0.360592*eavg1_mu*ferm2/ferm3 - 4.0*temp*ferm3_eta/ferm4_eta) + 0.31*ferm2_negeta*(eavg1_mu*ferm2/ferm3 - 4.0*temp*ferm3_negeta/ferm4_negeta)); 
  Real h2_mubar = 2.0*f1*lmu_y*(ferm2_eta*(0.31*eavg1_mu*ferm2/ferm3 - 4.0*temp*ferm3_eta/ferm4_eta) + 0.360592*ferm2_negeta*(eavg1_mu*ferm2/ferm3 - 4.0*temp*ferm3_negeta/ferm4_negeta)); 

  Real h2 = 1e12*(h2_nubar + h2_nu + h2_mu + h2_mubar)*(1.0-x)*inv_r2;
//neutrino nucleon scattering
   Real mn = 1.6726e-24;
   Real fact = (lnubar_y*eavg1_nubar*eavg1_nubar*(eavg1_nubar*ferm2/ferm3 - 6.0*temp*ferm5/ferm6) + lnu_y*eavg1_nu*eavg1_nu*(eavg1_nu*ferm2/ferm3 - 6.0*temp*ferm5/ferm6) + 4.0*lmu_y*eavg1_mu*eavg1_mu*(eavg1_mu*ferm2/ferm3 - 6.0*temp*ferm5/ferm6))/mn; //4.0 for mu, mubar, tau and taubar
//scattering on neutrons
   Real h3_n = 5.2225e4*(1.0-ye)*fact;
//scattering on protons
   Real h3_p = 4.3215e4*ye*fact;
   Real h3 = 6.2415e5*(h3_n + h3_p)*(1.0-x)*inv_r2;

//neutrino annihilation
  Real psi = (x*x+4.0*x+5.0)*std::pow(1.0-x,4.0);  
  Real h4 = (1.6e19*6.2415e5*1e32*psi/(rho*std::pow(r_0,4.0)))*(lnu_y*lnubar_y*(eavg1_nubar + eavg1_nu) + (6.0/7.0)*lmu_y*lmu_y*eavg1_mu);



  if(std::isnan(h1)) {
         h1 = 0.0; 
  }

  if(std::isnan(h2)) {
         h2 = 0.0; 
  } 

  if(std::isnan(h3)) {
         h3 = 0.0; 
  }

  if(std::isnan(h4)) {
         h4 = 0.0; 
  }

  Real heating = h1 + h2 + h3 + h4;


  Real ferm5_eta = fermi_approx(5,etaqw);
  Real ferm5_negeta = fermi_approx(5,-1.0*etaqw);
  Real cf = ye*(ferm5_eta/ferm5) + (1.0-ye)*(ferm5_negeta/ferm5);


  Real c1= 2.073*Na*cf*std::pow(temp,6.0); // cooling

//electron-positron annihilation
  Real c2 = (0.145*Na*1e8*std::pow(temp,9.0)/rho)*((ferm4_eta*ferm3_negeta + ferm4_negeta*ferm3_eta)/(2.0*ferm4*ferm3));

  if(std::isnan(c1)) {
         c1 = 0.0; 
  } 

  if(std::isnan(c2)) {
         c2 = 0.0; 
  } 


  Real cooling = c1 + c2; 

  if (return_index==0) {return heating - cooling;}
  else if (return_index==1) {return cooling/heating;}
  else {
    std::stringstream msg;
    msg << "Error in qdotQW" << std::endl
        << "return index different from 0 or 1 passed" << std::endl;
    ATHENA_ERROR(msg);
  }
}



inline Real zeroQ_temp(Real x, Real ye, Real time) {


  AthenaArray<Real> le;
  le.NewAthenaArray(12);
  lum_func(time,le);

  Real Coeff_nu = le(10);
  Real Coeff_nubar = le(11);


  return std::pow(1e12*(9.65/2.27)*((1.-ye)*Coeff_nu+ye*Coeff_nubar)*(1.-x)*inv_r2,
                  1./6.) / 8.6173e-11;
}



// exists_test1 from https://stackoverflow.com/a/12774387/2275975
inline bool exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

// Boundary Condition

void HSEInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                FaceField &b, Real time, Real dt, int is, int ie, int js,
                int je, int ks, int ke, int ngh);

void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);
                   
                   

//------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if(pin->GetInteger("mesh","x2rat")==-1) {
      std::cout<<"ENETR x2 space \n";
      EnrollUserMeshGenerator(X2DIR, grid_space_X2);
  }
  EnrollUserExplicitSourceFunction(UserSourceTerm);
  // set outer BC
  if (pin->GetString("mesh", "ox1_bc").compare("user") == 0) {
    if (Globals::my_rank == 0)
      printf("Using USER outer BC.\n");
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OutflowOuterX1);
  }
  // set inner BC
  int inner_BC_choice = pin->GetOrAddInteger("problem", "inner_BC_choice", 0);
  if (inner_BC_choice == 2) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HSEInnerX1);
  } 
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
        << "Invalid inner BC choice (" << inner_BC_choice << ")." << std::endl;
    ATHENA_ERROR(msg);
  }
 
  if (Globals::my_rank == 0) {printf("Using USER inner BC %d.\n", inner_BC_choice);}
  a_iso = pin->GetReal("hydro","iso_sound_speed");
  mu = pin->GetReal("problem","GM");
  Mpns = pin->GetOrAddReal("problem","Mpns",1.4);
  rho_0 = pin->GetReal("problem","rho_0");
  r_0 = pin->GetReal("mesh","x1min");
  inv_r2 = std::pow(r_0, -2.0);
  r0sq=r_0*r_0;
  Na = pin->GetReal("problem","Na");
  Ye=-100.0;
  ye_index=-1;
  if (pin->DoesParameterExist("hydro", "eos_ye")) {
    Ye = pin->GetReal("hydro", "eos_ye");
  }
  if (Globals::my_rank == 0) {std::cout<<"Ye   "<<Ye<<"\n";}
  if (pin->DoesParameterExist("hydro", "eos_ye_index")) {
    ye_index = pin->GetInteger("hydro", "eos_ye_index");
    if (ye_index < 0 || ye_index >= NSCALARS) {
      std::stringstream msg;
      msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
          << "hydro/eos_ye_index must be between 0 and NSCALARS (" << NSCALARS << ")."
          << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  if (Ye < 0 && ye_index < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
        << "either hydro/eos_ye or hydro/eos_ye_index must be specified."
        << std::endl;
    ATHENA_ERROR(msg);
  }

temp_index=-1;
if (pin->DoesParameterExist("hydro", "helm_temp_index")) {
      temp_index = pin->GetInteger("hydro", "helm_temp_index");
      if (Globals::my_rank == 0) {std::cout<<"Using Temp passive scalars in PGEN\n";}
      if (temp_index < 0 || temp_index >= NSCALARS) {
        std::stringstream msg;
        msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
            << "hydro/helm_temp_index must be between 0 and NSCALARS ("
            << NSCALARS << ")."
            << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  GuessTemp_Bc = 5.0;
  GuessYe_Bc = 0.05;
  fixed_tempBC = GuessTemp_Bc;
  Tc = pin->GetReal("problem","T_cutoff");
  B_0 = pin->GetReal("problem","B_0");
  B_0 = B_0/(std::pow(4.0*PI,0.5)); //convert to Lorentz-Heaviside units
  Omega_0 = pin->GetOrAddReal("problem","Omega_0",0.0);
  Omega_dot = pin->GetOrAddReal("problem","Omega_dot",0.0);
  alpha = pin->GetReal("problem","alpha");
  set_T_at_r0 = pin->GetOrAddBoolean("problem", "set_T_at_r0", false);
  rot_flag = pin->GetOrAddBoolean("problem", "rot_flag", false);
  rot_option = pin->GetInteger("problem", "rot_option");
  source_choice = pin->GetOrAddInteger("problem", "source_choice", 0);
  vr_choice = pin->GetOrAddInteger("problem", "vr_choice", 0);
  if (Globals::my_rank == 0) {std::cout<<"vr choice   "<<vr_choice<<"\n";}
  vth_choice = pin->GetOrAddInteger("problem", "vth_choice", 0);
  if (Globals::my_rank == 0) {std::cout<<"vth choice   "<<vth_choice<<"\n";}
  cs = pin->GetOrAddBoolean("problem", "cs", false);

  // initial limonosity/energies
  L_nu = pin->GetReal("problem","L_nu");
  L_mu = pin->GetReal("problem","L_mu");
  L_nubar = pin->GetReal("problem","L_nubar");
  eps_nu = pin->GetReal("problem","eps_nu");
  eps_mu = pin->GetReal("problem","eps_mu");
  eps_nubar = pin->GetReal("problem","eps_nubar");
  if (Globals::my_rank == 0) {std::cout<<"Lnu  "<<L_nu<<"  Lnubar   "<<L_nubar<<"   Lmu    "<<L_mu<<"  eps_nu   "<<eps_nu<<"  eps_nubar   "<<eps_nubar<<"  eps_mu   "<<eps_mu<<"\n";}

  Real ferm6 = 714.668;
  Real ferm5 = 118.266;
  Real ferm4 = 23.3309;
  Real ferm3 = 5.6822;
  Real ferm2 = 1.80309;
  
  Real ktemp_nubar = eps_nubar*ferm2/ferm3;
  Real ktemp_nu = eps_nu*ferm2/ferm3;
  Coeff_nu_0 = L_nu * SQR(ktemp_nu)*ferm5/ferm3;
  Coeff_nubar_0 = L_nubar * SQR(ktemp_nubar)*ferm5/ferm3;
  // finial limonosity/energies
  L_nu_f = pin->GetOrAddReal("problem","L_nu_f",L_nu);
  L_mu_f = pin->GetOrAddReal("problem","L_mu_f",L_mu);
  L_nubar_f = pin->GetOrAddReal("problem","L_nubar_f",L_nubar);
  eps_nu_f = pin->GetOrAddReal("problem","eps_nu_f",eps_nu);
  eps_mu_f = pin->GetOrAddReal("problem","eps_mu_f",eps_mu);
  eps_nubar_f = pin->GetOrAddReal("problem","eps_nubar_f",eps_nubar);

  Real ktemp_nubar_f = eps_nubar_f*ferm2/ferm3;
  Real ktemp_nu_f = eps_nu_f*ferm2/ferm3;

  Coeff_nu_f = L_nu_f * SQR(ktemp_nu_f)*ferm5/ferm3;
  dCoeff_nu = Coeff_nu_f - Coeff_nu_0;
  Coeff_nubar_f = L_nubar_f * SQR(ktemp_nubar_f)*ferm5/ferm3;
  dCoeff_nubar = Coeff_nubar_f - Coeff_nubar_0;
  Real inf = std::numeric_limits<Real>::infinity();
  t_L_0 = pin->GetOrAddReal("problem","t_L_0",inf);
  t_L_1 = pin->GetOrAddReal("problem","t_L_1",inf);
  if (t_L_1 < inf && t_L_1 <= t_L_0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
        << "l_transition_end <= l_transition_start" << std::endl;
    ATHENA_ERROR(msg);
  }

  //T_floor = pin->GetOrAddReal("hydro", "T_floor", float_eps);

  // Parse IC choice
  std::string file;
  bool has_file = pin->DoesParameterExist("problem", "file");
  if (has_file) {
    file = pin->GetString("problem", "file");
  }
  bool use_IC_specified = pin->DoesParameterExist("problem", "use_IC_file");
  use_IC_file = pin->GetOrAddBoolean("problem", "use_IC_file", has_file);

  if (use_IC_specified && use_IC_file) {
    if (!has_file) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
          << "No IC file specified in input file." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (!exists(file)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
          << "Specified IC file " << file << "does not exits." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  if (has_file) {
    if (!exists(file)) {
      use_IC_file = false;
      if (Globals::my_rank == 0) {
        std::cout << "Unable to locate IC file " << file << ", reverting to code IC."
                  << std::endl;
      }
    }
  }

  // Read ICs from data file
  if (use_IC_file) {
    rows = pin->GetInteger("problem", "rows");
    int cols = pin->GetInteger("problem", "cols");
    int col_rho = pin->GetInteger("problem", "col_rho");
    int col_v = pin->GetInteger("problem", "col_v");
    int col_T = pin->GetInteger("problem", "col_T");
    int col_ye = pin->GetInteger("problem", "col_ye");


    // Prepare arrays to hold profile
    AllocateRealUserMeshDataField(8);
    ruser_mesh_data[0].NewAthenaArray(rows);
    ruser_mesh_data[1].NewAthenaArray(rows);
    ruser_mesh_data[2].NewAthenaArray(rows);
    ruser_mesh_data[3].NewAthenaArray(rows);
    ruser_mesh_data[4].NewAthenaArray(1);
    ruser_mesh_data[5].NewAthenaArray(1);
    ruser_mesh_data[6].NewAthenaArray(1);
    ruser_mesh_data[7].NewAthenaArray(rows);
    AthenaArray<Real>& r_in{ruser_mesh_data[0]};
    AthenaArray<Real>& rho_in{ruser_mesh_data[1]};
    AthenaArray<Real>& v_in{ruser_mesh_data[2]};
    AthenaArray<Real>& T_in{ruser_mesh_data[3]};
    AthenaArray<Real>& Ye_in{ruser_mesh_data[7]};
    AthenaArray<Real>& Omega_in{ruser_mesh_data[5]};
    AthenaArray<Real>& rot_in{ruser_mesh_data[6]};
    Omega_in(0)=Omega_0;
    if(rot_flag){
      rot_in(0)=1;
    }
    else{
      rot_in(0)=0;
    }
    if (Globals::my_rank == 0)
      std::cout<< "Using IC file: " << file << "\n";
    std::string line;
    std::ifstream stream;

    stream.open(file);
    AthenaArray<Real> s_vals;
    s_vals.NewAthenaArray(cols);

    for (int n = 0; n < rows; ++n) {
      std::getline(stream, line);
      std::string word;
      std::stringstream iss(line);
      int m=0;
      while (iss >> word) {
        s_vals(m)=std::stof(word);
        m=m+1;
      }
      //std::cout<<line;
      r_in(n)=s_vals(0);
      rho_in(n) = s_vals(col_rho+1);
      v_in(n) = s_vals(col_v+1);
      T_in(n) = s_vals(col_T+1);
      Ye_in(n) = s_vals(col_ye+1);

      if (Globals::my_rank == 0) {
        //std::cout<<r_in(n)<<" ,"<<rho_in(n)<<" , "<<v_in(n)<<" , "<<T_in(n);
        //std::cout<<"\n";
      }
    }
  }


  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  if (lid == 0) {
    // Set up MPI communicator for inner r boundary conditions
    bool include_rank = false;
    int ntot = pmy_mesh->nblocal;//my_blocks.GetSize();
    for (int i = 0; i < ntot; ++i) {
      if (loc.lx1 == 0) {
        //std::cout << "Rank " << Globals::my_rank << " is in bc_group\n";
        include_rank = true;
        break;
      }
    }
    // Gather the inclusion information from all processes
    bool* include_ranks = new bool[Globals::nranks];
    for (int i = 0; i < Globals::nranks; ++i) {
      include_ranks[i] = false;
    }
    MPI_Allgather(&include_rank, 1, MPI_C_BOOL, include_ranks, 1, MPI_C_BOOL, MPI_COMM_WORLD);
    // Count the number of included ranks
    int n = 0;
    for (int i = 0; i < Globals::nranks; ++i) {
        if (include_ranks[i]) {
            ++n;
        }
    }
    // Create the ranks array for the new group
    int* ranks = new int[Globals::nranks];
    for (int i = 0; i < Globals::nranks; ++i) {
      ranks[i] = -1;
    }
    int index = 0;
    for (int i = 0; i < Globals::nranks; ++i) {
        if (include_ranks[i]) {
            ranks[index++] = i;
        }
    }
    n_bc_ranks = index;
    int expected = (pmy_mesh->mesh_size.nx2 / block_size.nx2) * (pmy_mesh->mesh_size.nx3 / block_size.nx3);
    if (!Globals::my_rank) {
      //std::cout << "Total ranks in BC comm: " << n_bc_ranks << ", " << n << ", " << expected << "\n";
    }
    if (n_bc_ranks != n || n != expected) {
      std::stringstream msg;
      msg << "";
      if (Globals::my_rank == 0) {
        msg << "### FATAL ERROR in InitUserMeshBlockData" << std::endl
            << " BC block count mismatch: " << n << ", " << n_bc_ranks << ", "
            << expected << std::endl;
      }
      ATHENA_ERROR(msg);
    }
    delete[] include_ranks;
    if (n_bc_ranks) {
      // Create a group from the world communicator
      MPI_Group world_group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);

      // Create a new group based on the inclusion list
      MPI_Group new_group;
      MPI_Group_incl(world_group, n_bc_ranks, ranks, &new_group);

      // Create a new communicator based on the new group
      MPI_Comm_create(MPI_COMM_WORLD, new_group, &bc_comm);
    }
    delete[] ranks;
  }
  AllocateIntUserMeshBlockDataField(1);
  iuser_meshblock_data[0].NewAthenaArray(1);
  int& index_compute{iuser_meshblock_data[0](0)};
  index_compute = -1;
  Real r_compute = pin->GetOrAddReal("problem", "r_compute", 0.0);
  if (block_size.x1min < r_compute && r_compute <= block_size.x1max) {
    for (int i=is; i<=ie; i++) {
      // if outer radial face of cell > r_compute
      if (pcoord->x1f(i+1) > r_compute) {
        index_compute = i;
        break;
      }
    }
  }
  if (GENERAL_EOS) {
         ;
  }
  else {
    p_0=a_iso*a_iso*rho_0;
  }
      
  int i=0;
  if (COMP_DT) {

      IDT1 = i;
      IDT2 = i + 1;
      IDT3 = i + 2;
      IDT4 = i + 3;
      IDT5 = i + 4;
      i += 5;
      i=i+4;
  }
  if (COMP_DIV_B) {
      i+=1;
  }
  if (COMP_SRC && rot_flag) {
//some correctionns required for uov saving in rot_frame
    ISRC1 = i;
    ISRC2 = i + 1;
    ISRC3 = i + 2;
    ISRC4 = i + 3;
    ISRC5 = i + 4;
    ISRC6 = i + 5;
    ISRC7 = i + 6;
    ISRC8 = i + 7;
    ISRC9 = i + 8;
    i += 9;
  }
  else {
    i+=1;
  }

  if (i) AllocateUserOutputVariables(i);

  SetUserOutputVariableName(IDT1, "dt1"); //Temp
  SetUserOutputVariableName(IDT2, "dt2"); //qdot
  SetUserOutputVariableName(IDT3, "dt3"); //cs
  SetUserOutputVariableName(IDT4, "l_nu");
  SetUserOutputVariableName(IDT5, "l_nubar");
  SetUserOutputVariableName(5, "l_mu");
  SetUserOutputVariableName(6, "eps_nu");
  SetUserOutputVariableName(7, "eps_nubar");
  SetUserOutputVariableName(8, "eps_mu");
  if (COMP_DIV_B) {
    SetUserOutputVariableName(10, "divB");
  }    
  if (COMP_SRC && rot_flag) {
    SetUserOutputVariableName(ISRC1, "cor1");
    SetUserOutputVariableName(ISRC2, "cor2");
    SetUserOutputVariableName(ISRC3, "cor3");
    SetUserOutputVariableName(ISRC4, "cen1");
    SetUserOutputVariableName(ISRC5, "cen2");
    SetUserOutputVariableName(ISRC6, "cen3");
    SetUserOutputVariableName(ISRC7, "eng1");
    SetUserOutputVariableName(ISRC8, "eng2");
    SetUserOutputVariableName(ISRC9, "Omega");
  }
  else {
    SetUserOutputVariableName(9, "Omega");
  }
  
}

static Real ConservedRotationProfile(const Real x1, const Real x2) {
  if(rot_flag){
    return Omega_0 * std::sin(x2) * x1 * (std::pow((r_0/x1),2.0) - 1.0);
  }
  else {
    return Omega_0 * std::sin(x2) * x1 * (std::pow((r_0/x1),2.0));
  }
}




// Finds the QW EoS electron chemical potential.
Real QWEta(Real rho, Real T, Real Ye) {
  // Returns eta = mu_e / kbT
  // Physical constants
  Real third         = 1.0 / 3.0;
  Real c             = 2.99792458e10;              // Speed of light in cm/s
  Real k             = 1.380649e-16;               // Boltzmann constant in erg/K
  Real mn            = 1.6726e-24;                 // Baryon mass in g
  Real hbar          = 6.62607015e-27/(2.0*PI);    // Reduced Planck's constant in erg s
  Real c3            = std::pow(k/(hbar*c),3.0);
  Real eta_den_const = std::pow(6.0,2.0*third);
  Real root3         = std::sqrt(3.0);
  Real eta_den_A     = std::pow(2.0,third) / eta_den_const;
  Real eta_den_B     = 2.0*std::pow(3.0,third)/eta_den_const;

  Real vol       = mn/rho;
  Real T3        = T*T*T;
  Real T4        = T*T3;
  Real a         = c3*std::pow(T,3.0)*vol*(PI/3.0);
  Real a2        = SQR(a);
  Real a4        = SQR(a2);
  Real a6        = a2 * a4;
  Real y2        = SQR(Ye);
  Real b         = std::sqrt(4.0*a6+27.0*a4*y2);
  Real term      = std::pow(9.0*a2*Ye+root3*b, third);
  Real eta_by_pi = (eta_den_A)*term/a - (eta_den_B)*a/term; // actually eta/pi
  Real eta       = eta_by_pi * PI;
  if ((std::isnan(T)) || (std::isnan(rho)) || (std::isnan(eta)) || (std::isnan(Ye))) {
      
      std::stringstream msg;
      msg <<"Nan in PGEN, QWEta \n";
      ATHENA_ERROR(msg);

  }
  
  if(eta<0.0) {
      std::stringstream msg;
      msg <<"Negative ETA being returned in PGEN QWEta \n";
      ATHENA_ERROR(msg);
  }
  return eta;

}


// Ye source function for scalar Ye
Real YeSource(Real temp, Real ye, Real x, Real time, Real Rho, AthenaArray<Real> &le, int return_index) {
  //return_index 0 for main and 1 for BCs

  // Returns Sources for Yedot=vdYe/dr in s^-1 units
  // Temperature argument is in MeV
  // Define constants
  Real alpha    = 1.254;                                // Coupling Coefficient
  Real G_F      = 1.16637e-11; // Fermi Constant in MeV^-2
  Real Delta    = 1.2935;                               // neutron-proton mass difference in MeV
  Real hbar     = 6.582119569e-22;                     // Planck's constant in MeV s
  Real c        = 2.99792458e10;                       // Light speed in cm/s
  Real erg2MeV  = 6.24151e5;                           // convert erg to MeV

//Here epsnu_y and epsnubar_y are square roots of the energy term in charged current heating rate
  Real lnu_y = le(0);
  Real lnubar_y = le(1);
  Real lmu_y = le(2);
  Real epsnu_y = le(3);
  Real epsnubar_y = le(4);
  Real eavg1_mu = le(5);
  Real eavg1_nu = le(6) ;
  Real eavg2_nu = le(7);
  Real eavg1_nubar = le(8);
  Real eavg2_nubar = le(9);
  Real Coeff_nu = le(10);
  Real Coeff_nubar = le(11);

  Real ferm6 = 714.668;
  Real ferm5 = 118.266;
  Real ferm4 = 23.3309;
  Real ferm3 = 5.6822;
  Real ferm2 = 1.80309;
  
  Real lambda_nue_n  = std::pow(hbar*c,2.0)*((1.0+3.0*alpha*alpha)/(2.0*PI*PI))*G_F*G_F*lnu_y*(1e51*erg2MeV/(r_0*r_0))*(eavg2_nu/eavg1_nu+2.0*Delta+Delta*Delta/eavg1_nu)*(1.0-x); // s^-1 units
  Real lambda_nuebar_p  = std::pow(hbar*c,2.0)*((1.0+3.0*alpha*alpha)/(2.0*PI*PI))*G_F*G_F*lnubar_y*(1e51*erg2MeV/(r_0*r_0))*(eavg2_nubar/eavg1_nubar-2.0*Delta+Delta*Delta/eavg1_nubar)*(1.0-x); // s^-1 units

// Calculate Eta, send in T in K
  Real Eta = QWEta(Rho,temp/8.6173e-11,ye);

  // Reverse reaction rates
  Real f4_eta       = fermi_approx(4,Eta);
  Real f4_negeta    = fermi_approx(4,-1.0*Eta);
  Real cc = 0.448*std::pow(temp,5.0)/ferm4;
  Real lambda_ele_p = cc*f4_eta; // s^-1 units (e- p)
  Real lambda_pos_n = cc*f4_negeta;          // s^-1 units (e+ n)



  if(std::isnan(lambda_nue_n)) {
         lambda_nue_n = 0.0; 
  }



  if(std::isnan(lambda_nuebar_p)) {
         lambda_nuebar_p = 0.0; 
  }


  if(std::isnan(lambda_ele_p)) {
         lambda_ele_p = 0.0; 
  }

  if(std::isnan(lambda_pos_n)) {
         lambda_pos_n = 0.0; 
  } 


  // Source term
  //std::cout<<"ENTER YeSOurce   "<<L_nu<<"    "<<eps_nu<<"    "<<lambda_nue_n<<"   "<<lambda_pos_n<<"    "<<lambda_ele_p<<"   "<<lambda_nuebar_p<<"    "<<temp/8.6173e-11<<"    "<<ye<<"    "<<Rho<<"\n";
  if(return_index==0) {return lambda_nue_n+lambda_pos_n-(lambda_nue_n+lambda_pos_n+lambda_nuebar_p+lambda_ele_p)*ye;}
  else if(return_index==1) {return (lambda_nue_n+lambda_pos_n+lambda_nuebar_p+lambda_ele_p)*ye/(lambda_nue_n+lambda_pos_n);}
  
  else {
    std::stringstream msg;
    msg << "Error in YeSource" << std::endl
        << "return index different from 0 or 1 passed" << std::endl;
    ATHENA_ERROR(msg);
  }

}







void BC_T0(Real Ye, Real guessT, Real time, Real r, Real theta, AthenaArray<Real> &out) {

  AthenaArray<Real> le;
  le.NewAthenaArray(12);
  lum_func(time,le);
  //std::cout<<"Lums   "<<le(0)<<"     "<<le(1)<<"\n";
  Real kbol = 1.380649e-16;
  Real mn = 1.6726e-24;
  Real kmn = kbol/(mn*8.6173e-11);
  Real omega_t = Omega_0;

  Real r0 = r_0;

  Real e1 = (r0 - r) * mu / (r * r0);
  Real e2 = 0.5*omega_t*omega_t;


  Real e3 = r*r*std::sin(theta)*std::sin(theta);
  

  Real dpdd_0 = kmn*guessT;
  Real rho = rho_0*std::exp(e1/dpdd_0) * std::exp((e2/dpdd_0)*(e3-r0sq));
  Real eta = QWEta(rho,guessT/8.6173e-11,Ye);
  Real ht = 1e-4;
  Real hy = 1e-4;
  Real T2 = guessT + ht*guessT;

  Real dpdd_0_T2 = kmn*T2;
  Real rho_T2 = rho_0*std::exp(e1/dpdd_0_T2) * std::exp((e2/dpdd_0_T2)*(e3-r0sq));
  Real eta_T2 = QWEta(rho_T2,T2/8.6173e-11,Ye); 

 
  Real Ye2 = Ye + hy*Ye;
  Real eta_Y2 = QWEta(rho,guessT/8.6173e-11,Ye2); 

  Real heat = qdotQW(guessT, rho,0.0, Ye, time, eta, le,1);
  Real qratm = heat - 1.0;
  Real qrat_dT = (qdotQW(T2,rho_T2, 0.0, Ye, time, eta_T2, le,1) - 1.0 -qratm)/(ht*guessT);
  Real qrat_dY = (qdotQW(guessT,rho, 0.0, Ye2, time, eta_Y2, le,1) - 1.0 -qratm)/(hy*Ye);

  Real Yratm = YeSource(guessT, Ye, 0.0, time, rho, le,1) - 1.0; 
  Real Yrat_dT = (YeSource(T2, Ye, 0.0, time, rho_T2, le,1) - 1.0 - Yratm)/(ht*guessT); 
  Real Yrat_dY = (YeSource(guessT, Ye2, 0.0, time, rho, le,1) - 1.0 - Yratm)/(hy*Ye); 



  int nlim = 100000;
  int iter=0;
  Real prec_bc=1e-4;
  Real step;
  //std::cout<<"In BC_T0   "<<guessT<<"     "<<Ye<<"     "<<rho<<"     "<<eta<<"     "<<"      "<<qrat_dT<<"      "<<qrat_dY<<"      "<<Yrat_dT<<"      "<<Yrat_dY<<"      "<<qratm<<"     "<<Yratm<<"\n";
  while ((std::abs(qratm)>=prec_bc) || (std::abs(Yratm)>=prec_bc)) {  
    iter = iter+1;

 //inverse of jacobian

    Real a = qrat_dT;
    Real b = qrat_dY;
    Real c = Yrat_dT;
    Real d = Yrat_dY;

    Real det = a*d-b*c;
    Real ai = d/det;
    Real bi = -1.0*b/det;
    Real ci = -1.0*c/det;
    Real di = a/det;

    Real iT = ai*qratm + bi*Yratm;
    Real iY = ci*qratm + di*Yratm;
    
    step = 1.0;
    guessT = guessT - step*iT;
    Ye = Ye - step*iY;

    dpdd_0 = kmn*guessT;
    rho = rho_0*std::exp(e1/dpdd_0) * std::exp((e2/dpdd_0)*(e3-r0sq));
    eta = QWEta(rho,guessT/8.6173e-11,Ye);
    T2 = guessT + ht*guessT;

    dpdd_0_T2 = kmn*T2;
    rho_T2 = rho_0*std::exp(e1/dpdd_0_T2) * std::exp((e2/dpdd_0_T2)*(e3-r0sq));
    eta_T2 = QWEta(rho_T2,T2/8.6173e-11,Ye); 

 
    Ye2 = Ye + hy*Ye;
    eta_Y2 = QWEta(rho,guessT/8.6173e-11,Ye2); 

    heat = qdotQW(guessT, rho,0.0, Ye, time, eta, le,1);
    qratm = heat - 1.0;
    qrat_dT = (qdotQW(T2,rho_T2, 0.0, Ye, time, eta_T2, le,1) - 1.0 -qratm)/(ht*guessT);
    qrat_dY = (qdotQW(guessT,rho, 0.0, Ye2, time, eta_Y2, le,1) - 1.0 -qratm)/(hy*Ye);

    Yratm = YeSource(guessT, Ye, 0.0, time, rho, le,1) - 1.0; 
    Yrat_dT = (YeSource(T2, Ye, 0.0, time, rho_T2, le,1) - 1.0 - Yratm)/(ht*guessT); 
    Yrat_dY = (YeSource(guessT, Ye2, 0.0, time, rho, le,1) - 1.0 - Yratm)/(hy*Ye); 

  //std::cout<<"In BC_T0   "<<guessT<<"     "<<Ye<<"     "<<rho<<"     "<<eta<<"     "<<"      "<<qrat_dT<<"      "<<qrat_dY<<"      "<<Yrat_dT<<"      "<<Yrat_dY<<"     "<<qratm<<"     "<<Yratm<<"\n";
    if ((std::isnan(guessT)) || (std::isnan(rho)) || (std::isnan(Ye))) {
      
      std::stringstream msg;
      msg <<"Nan in PGEN, BC_T0 \n";
      ATHENA_ERROR(msg);
   }
    if (iter>nlim) {
      
      std::stringstream msg;
      msg <<"Nlim exceeded in BC_T0 \n";
      ATHENA_ERROR(msg);
   }

 }

  if ((guessT < 0.0) || (Ye < 0.0)) {
      
      std::stringstream msg;
      msg <<"Negative Temp or Ye being returned in BC_T0 \n";
      ATHENA_ERROR(msg);
   }
   

  if ((guessT < 0.1) || (guessT > 10.0) || (Ye > 0.2)) {
      std::cout<<guessT<<"        "<<Ye<<"\n";
      std::stringstream msg;
      msg <<"Very small or very large T, Ye being returned in BC_T0 \n";
      ATHENA_ERROR(msg);
   }
  //std::cout<<"In BC_T0   "<<guessT<<"     "<<Ye<<"     "<<rho<<"     "<<eta<<"     "<<"      "<<qrat_dT<<"      "<<qrat_dY<<"      "<<Yrat_dT<<"      "<<Yrat_dY<<"     "<<qratm<<"     "<<Yratm<<"\n";
  //std::cout<<"In BC_T0   "<<guessT<<"     "<<Ye<<"     "<<iter<<"\n";
  //std::cout<<"Returning In BC_T0   "<<r<<"     "<<guessT<<"     "<<Ye<<"     "<<rho_0<<"     "<<rho<<"     "<<"      "<<le(0)<<"      "<<le(1)<<"      "<<le(2)<<"      "<<le(5)<<"     "<<le(6)<<"     "<<le(8)<<"\n";
  out(0) = guessT;
  out(1) = Ye;

} 



void BC_T0_gh2(Real Ye, Real guessT, Real time, Real r, Real theta, AthenaArray<Real> &out) {

  AthenaArray<Real> le;
  le.NewAthenaArray(12);
  lum_func(time,le);
  //std::cout<<"Lums   "<<le(0)<<"     "<<le(1)<<"\n";
  Real kbol = 1.380649e-16;
  Real mn = 1.6726e-24;
  Real kmn = kbol/(mn*8.6173e-11);
  Real omega_t = Omega_0;

  Real r0 = r_0;

  Real e1 = (r0 - r) * mu / (r * r0);
  Real e2 = 0.5*omega_t*omega_t;

  Real e3 = r*r*std::sin(theta)*std::sin(theta);
  

  Real dpdd_0 = kmn*guessT;
  Real rho = rho_0*std::exp(e1/dpdd_0) * std::exp((e2/dpdd_0)*(e3-r0sq));
  Real eta = QWEta(rho,guessT/8.6173e-11,Ye);
  Real hy = 1e-4;
  
 
  Real Ye2 = Ye + hy*Ye;
  Real eta_Y2 = QWEta(rho,guessT/8.6173e-11,Ye2); 

  Real heat = qdotQW(guessT, rho,0.0, Ye, time, eta, le,1);
  Real qratm = heat - 1.0;

  Real qrat_dY = (qdotQW(guessT,rho, 0.0, Ye2, time, eta_Y2, le,1) - 1.0 -qratm)/(hy*Ye);

  Real Yratm = YeSource(guessT, Ye, 0.0, time, rho, le,1) - 1.0; 

 


  int nlim = 100000;
  int iter=0;
  Real prec_bc=1e-4;
  Real step;
  //std::cout<<"In BC_T0   "<<guessT<<"     "<<Ye<<"     "<<rho<<"     "<<eta<<"     "<<"      "<<qrat_dT<<"      "<<qrat_dY<<"      "<<Yrat_dT<<"      "<<Yrat_dY<<"      "<<qratm<<"     "<<Yratm<<"\n";
   
  while ((std::abs(qratm)>=prec_bc)) {  
    iter = iter+1;


    Real iY = qratm/qrat_dY;
    step = 1.0;

    Ye = Ye - step*iY;

    eta = QWEta(rho,guessT/8.6173e-11,Ye);

 
    Ye2 = Ye + hy*Ye;
    eta_Y2 = QWEta(rho,guessT/8.6173e-11,Ye2); 

    heat = qdotQW(guessT, rho,0.0, Ye, time, eta, le,1);
    qratm = heat - 1.0;

    qrat_dY = (qdotQW(guessT,rho, 0.0, Ye2, time, eta_Y2, le,1) - 1.0 -qratm)/(hy*Ye);

    Yratm = YeSource(guessT, Ye, 0.0, time, rho, le,1) - 1.0; 

  //std::cout<<"In BC_T0   "<<guessT<<"     "<<Ye<<"     "<<rho<<"     "<<eta<<"     "<<"      "<<qrat_dT<<"      "<<qrat_dY<<"      "<<Yrat_dT<<"      "<<Yrat_dY<<"     "<<qratm<<"     "<<Yratm<<"\n";
    if ((std::isnan(guessT)) || (std::isnan(rho)) || (std::isnan(Ye))) {
      
      std::stringstream msg;
      msg <<"Nan in PGEN, BC_T0 gh2 \n";
      ATHENA_ERROR(msg);
   }
    if (iter>nlim) {
      
      std::stringstream msg;
      msg <<"Nlim exceeded in BC_T0 gh2\n";
      ATHENA_ERROR(msg);
   }

 }

  if ((guessT < 0.0) || (Ye < 0.0)) {
      
      std::stringstream msg;
      msg <<"Negative Temp or Ye being returned in BC_T0 gh2\n";
      ATHENA_ERROR(msg);
   }
   
   //if (std::abs(Yratm)>0.25) {
   
      //std::stringstream msg;
      //msg <<"Yedot precision greater than 25% in BC_T0 gh2\n";
      //ATHENA_ERROR(msg);
   //}
   
   if(Ye>0.2) {
      std::cout<<guessT<<"        "<<Ye<<"\n";
      std::stringstream msg;
      msg <<"Very large Ye being returned in BC_T0 gh2 \n";
      ATHENA_ERROR(msg);
   }
  //std::cout<<"In BC_T0  gh2 "<<guessT<<"     "<<Ye<<"     "<<rho<<"     "<<eta<<"     "<<"     "<<qratm<<"     "<<Yratm<<"\n";
  //std::cout<<"In BC_T0   "<<guessT<<"     "<<Ye<<"     "<<iter<<"\n";
  //std::cout<<"Returning In BC_T0  gh2 "<<r<<"     "<<guessT<<"     "<<Ye<<"     "<<rho_0<<"     "<<rho<<"     "<<"      "<<le(0)<<"      "<<le(1)<<"      "<<le(2)<<"      "<<le(5)<<"     "<<le(6)<<"     "<<le(8)<<"\n";
  out(0) = guessT;
  out(1) = Ye;

}



//------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Parker wind

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // if root processor and zeroth block
  if ((Globals::my_rank == 0) && (lid == 0)) {
    std::cout<<"ENTER PGEN \n";
    std::cout<<NGHOST<<"   "<<is<<"   "<<pcoord->x1f(is)<<"\n";;
  }

  if (use_IC_file) {
    // define references for MeshBlock::ProblemGenerator
    AthenaArray<Real>& r_in{pmy_mesh->ruser_mesh_data[0]};
    AthenaArray<Real>& rho_in{pmy_mesh->ruser_mesh_data[1]};
    AthenaArray<Real>& v_in{pmy_mesh->ruser_mesh_data[2]};
    AthenaArray<Real>& T_in{pmy_mesh->ruser_mesh_data[3]};
    AthenaArray<Real>& Ye_in{pmy_mesh->ruser_mesh_data[7]};


    for (int k=ks; k<=ke; k++) {

      // Real phi = pcoord->x3v(k);
      for (int j=js; j<=je; j++) {
	
        Real theta = pcoord->x2v(j);
	Real sinx2l = std::sin(pcoord->x2f(j));
	Real sinx2u = std::sin(pcoord->x2f(j+1));
	Real cosx2l = std::cos(pcoord->x2f(j));
	Real cosx2u = std::cos(pcoord->x2f(j+1));
	Real x2l = pcoord->x2f(j);
	Real x2u = pcoord->x2f(j+1);
	
	  
        for (int i=is; i<=ie; i++) {
          Real r = pcoord->x1v(i);
	      Real x1l = pcoord->x1f(i);
	      Real x1u = pcoord->x1f(i+1);
          Real rho, v, temp, vth, ye;

          int index=0;
          Real min=1e15;
          Real diff;

          for (int f=0; f<rows; f++) {
            diff=r-r_in(f);
            if(diff>=0.0) {
              if(diff<min) {
                min=diff;
                index=f;
              }
            }
          }
	  
          //linear interpolation when r values in ICs and simulation are different
	  

	    
          rho=rho_in(index)+(r-r_in(index))*(rho_in(index+1)-rho_in(index))
                    /(r_in(index+1)-r_in(index));
	   
          v=v_in(index)+(r-r_in(index))*(v_in(index+1)-v_in(index))
                  /(r_in(index+1)-r_in(index));
          

          temp=T_in(index)+(r-r_in(index))*(T_in(index+1)-T_in(index))
               /(r_in(index+1)-r_in(index));
          if (NSCALARS > 0 && ye_index >= 0) {
            ye=Ye_in(index)+(r-r_in(index))*(Ye_in(index+1)-Ye_in(index))
               /(r_in(index+1)-r_in(index));
          }
          else {
            ye=Ye;
          }
          //std::cout<<"ye in pgen func "<<ye<<"\n";
	  
          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = v * rho;
          phydro->u(IM2,k,j,i) = 0.0;
          if (ye_index>=0) {
             pscalars->s(ye_index,k,j,i) = ye * rho;

             pscalars->r(ye_index,k,j,i) = ye;

          }
          else {
             pscalars->s(0,k,j,i) = ye * rho;
             pscalars->r(0,k,j,i) = ye;
          }

          if (temp_index>=0) {
             pscalars->s(temp_index,k,j,i) = temp * rho;

             pscalars->r(temp_index,k,j,i) = temp;

          }
	  phydro->u(IM3,k,j,i) = ConservedRotationProfile(r,theta)*phydro->u(IDN,k,j,i);

          if (NON_BAROTROPIC_EOS) {
            if (GENERAL_EOS) {
            Real s_cell[NSCALARS];
            for (int n=0; n<NSCALARS; ++n) {
              s_cell[n] = pscalars->s(n,k,j,i);
            }
            Real r_cell[NSCALARS];
            for (int n=0; n<NSCALARS; ++n) {
              r_cell[n] = pscalars->r(n,k,j,i);
            }
            //std::cout<<"r_cell  "<<r_cell[0]<<"\n";
            Real pressure = peos->PresFromRhoT(rho, temp,r_cell);
            phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(rho, pressure,r_cell);
            }
            phydro->u(IEN,k,j,i)+=0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       +SQR(phydro->u(IM3,k,j,i)))/rho;
          }
        }
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    // if root processor and zeroth block
    if ((Globals::my_rank == 0) && (lid == 0)) {
      std::cout<<"YES ENTER B field\n";
    }
    AthenaArray<Real> a1,a2,a3;
    int nx1 = (ie-is)+1 + 2*(NGHOST);
    int nx2 = (je-js)+1 + 2*(NGHOST);
    int nx3 = (ke-ks)+1 + 2*(NGHOST);
    a1.NewAthenaArray(nx3,nx2,nx1);
    a2.NewAthenaArray(nx3,nx2,nx1);
    a3.NewAthenaArray(nx3,nx2,nx1);

    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
        }
      }
    }

    // Initialize interface fields
    AthenaArray<Real> area,len,len_p1;
    area.NewAthenaArray(nx1);
    len.NewAthenaArray(nx1);
    len_p1.NewAthenaArray(nx1);

    // for 1,2,3-D
    int jl=js; int ju=je+1;
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary

      if ((pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) ||
        (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge))
        jl=js+1;
      if ((pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ||
        (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
        ju=je;
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face2Area(k,j,is,ie,area);
        pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = -1.0*(len(i+1)*a3(k,j,i+1) - len(i)*a3(k,j,i))/area(i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->Face3Area(k,j,is,ie,area);
        pcoord->Edge2Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = (len(i+1)*a2(k,j,i+1) - len(i)*a2(k,j,i))/area(i);
        }
      }
    }

    // for 2D and 3D
    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge3Length(k,j  ,is,ie+1,len);
          pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = (len_p1(i)*a3(k,j+1,i) - len(i)*a3(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face3Area(k,j,is,ie,area);
          pcoord->Edge1Length(k,j  ,is,ie,len);
          pcoord->Edge1Length(k,j+1,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) -= (len_p1(i)*a1(k,j+1,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
    // for 3D only
    if (block_size.nx3 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge2Length(k  ,j,is,ie+1,len);
          pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) -= (len_p1(i)*a2(k+1,j,i) - len(i)*a2(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke; ++k) {
        // reset loop limits for polar boundary
        int jl=js; int ju=je+1;
        if ((pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) ||
          (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge))
          jl=js+1;
        if ((pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ||
          (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
          ju=je;
        for (int j=jl; j<=ju; ++j) {
          pcoord->Face2Area(k,j,is,ie,area);
          pcoord->Edge1Length(k  ,j,is,ie,len);
          pcoord->Edge1Length(k+1,j,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) += (len_p1(i)*a1(k+1,j,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
    // Calculate cell-centered magnetic field
    AthenaArray<Real> bb;
    bb.NewAthenaArray(3, ke+1, je+1, ie+NGHOST+1);
    pfield->CalculateCellCenteredField(pfield->b, bb, pcoord, is-NGHOST, ie+NGHOST, js,
                                       je, ks, ke);

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real& bcc1 = bb(IB1,k,j,i);
          Real& bcc2 = bb(IB2,k,j,i);
          Real& bcc3 = bb(IB3,k,j,i);

          phydro->u(IEN,k,j,i) += 0.5*(SQR(bcc1)+SQR(bcc2)+SQR(bcc3));
         }
      }
    }
    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
    bb.DeleteAthenaArray();
  } // end if MAGNETIC_FIELDS_ENABLED
  // get energy floor
  //energy_floor = pin->GetReal("hydro", "efloor");

  // if root processor and last block
  if ((Globals::my_rank == 0) && (lid == pmy_mesh->nblocal - 1)) {
    std::cout<<"EXIT PGEN \n";
  }
}


// Source Terms
void heat_cool(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
               const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar) { 
  //std::cout<<"ENTER HEAT COOL \n";
  AthenaArray<Real> le;
  le.NewAthenaArray(12);
  lum_func(time,le);

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real r = pmb->pcoord->x1v(i);
        Real p = prim(IPR,k,j,i);
        Real rho = prim(IDN,k,j,i);
        Real r_cell[NSCALARS];
        for (int n=0; n<NSCALARS; ++n) {
           r_cell[n] = pmb->pscalars->r(n,k,j,i);
        }
        //std::cout<<"In heat_cool \n";
        
        Real temp = pmb->peos->TFromRhoP(rho, p, r_cell);
        
        Real x = std::sqrt(1.0-(r_0*r_0)/(r*r));
        if (NSCALARS > 0 && ye_index >= 0) {
           Ye = prim_scalar(ye_index,k,j,i);
        }
        Real Eta = QWEta(rho,temp,Ye);
        temp = temp * 8.6173e-11;  //T in MeV

        Real qdot = qdotQW(temp, rho, x, Ye, time, Eta,le,0); //MeV s^{-1} g^{-1} units
        //if(r<2e6) {std::cout<<"Qdot    "<<qdot*1.6022e-6/1e21<<"\n";}
       
        Real de = dt*prim(IDN,k,j,i) * qdot * 1.6022e-6; 
        cons(IEN,k,j,i) += de;
        if (NSCALARS > 0 && ye_index >= 0) {
           Real Ye_i = prim_scalar(ye_index,k,j,i);
	       Real dYe = dt*YeSource(temp,Ye_i,x,time,rho,le,0);
           //std::cout<<"dYe in heat_cool   "<<dYe<<"   "<<dt<<"\n";
           cons_scalar(ye_index,k,j,i) += dYe*prim(IDN,k,j,i);
        }

      }
    }
  }
  return;
}


//Rot source term

void RotFrame(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
               const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar) {
  //std::cout<<"ENTER ROT \n";
  // a_cor = -2 * Omega x v
  Real x1,x2,x3;
  Real f = dt * Omega_0;
  Real cen_dt = f * Omega_0;
  Real ff = f*f;
  AthenaArray<Real> &x1flux=pmb->phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[X2DIR];
  // AthenaArray<Real> &x3flux=pmb->phydro->flux[X3DIR];
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    Real dx3 = pmb->pcoord->dx3f(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      Real sinx2 = std::sin(x2);
      Real cosx2 = std::cos(x2);
      Real sin2x2 = std::sin(2 * x2);
      Real cos2x2 = std::cos(2 * x2);
      Real sinx2l = std::sin(pmb->pcoord->x2f(j));
      Real sinx2r = std::sin(pmb->pcoord->x2f(j+1));
      Real cosx2l = std::cos(pmb->pcoord->x2f(j));
      Real cosx2r = std::cos(pmb->pcoord->x2f(j+1));
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        x1 = pmb->pcoord->x1v(i);
        Real vol = pmb->pcoord->GetCellVolume(k,j,i);

        // Coriolis
        Real mom1 = cons(IM1,k,j,i);
        Real mom2 = cons(IM2,k,j,i);
        Real mom3 = cons(IM3,k,j,i);
        cons(IM1,k,j,i) += ff*(mom1*cos2x2 - mom2*sin2x2) + 2*f*mom3*sinx2;
        cons(IM1,k,j,i) /= (1+ff);
        cons(IM2,k,j,i) += -ff*(mom1*sin2x2 + mom2*cos2x2) + 2*f*mom3*cosx2;
        cons(IM2,k,j,i) /= (1+ff);
        cons(IM3,k,j,i) -= ff*mom3 + 2*f*(mom1*sinx2 + mom2*cosx2);
        cons(IM3,k,j,i) /= (1+ff);

        // Centrifugal
        Real dr4 = std::pow(pmb->pcoord->x1f(i+1),4) - std::pow(pmb->pcoord->x1f(i),4);
        Real cen_0 = 1./12. * dx3 / vol * cen_dt * dr4;
        Real src1 = .25*cen_0*(-9*(cosx2r-cosx2l)
                               -3*(cosx2r*SQR(sinx2r)-cosx2l*SQR(sinx2l))
                               +std::pow(cosx2r,3) -std::pow(cosx2l,3));
        cons(IM1,k,j,i) += src1 * prim(IDN,k,j,i);
        Real src2 = cen_0*(std::pow(sinx2r,3)-std::pow(sinx2l,3));
        cons(IM2,k,j,i) += src2 * prim(IDN,k,j,i);

        // Energy
        if (NON_BAROTROPIC_EOS) {
          // dE = a . v * dt
          cons(IEN,k,j,i) += src1 * 0.5*(x1flux(IDN,k,j,i) + x1flux(IDN,k,j,i+1));
          if (pmb->block_size.nx2 > 1) {
            cons(IEN,k,j,i) += src2 * 0.5*(x2flux(IDN,k,j,i) + x2flux(IDN,k,j+1,i));
          } else {
            cons(IEN,k,j,i) += src2 * prim(IDN,k,j,i) * prim(IVX,k,j,i);
          }
          // (src3=0)     += src3 * 0.5*(x3flux(IDN,k,j,i) + x3flux(IDN,k+1,j,i));
        }
      }
    }
  }
  return;
}

//Zhu rot source terms
void RotFrame2(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim,
              const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc,
              AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  // a_cor = -2 * Omega x v
  Real x1,x2,x3,sinx2;
  Real sm,sp,cm,cp,rm,rp;
  Real src1_i, src1_j;
  Real flux_kp,flux_km,flux_jp,flux_jm,flux_ip,flux_im;
  Real cor_[3];
  //~ AthenaArray<Real> &x1flux=pmb->phydro->flux[X1DIR];
  //~ AthenaArray<Real> &x2flux=pmb->phydro->flux[X2DIR];
  //~ AthenaArray<Real> &x3flux=pmb->phydro->flux[X3DIR];
  Real omega_t = Omega_0;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      if (pmb->block_size.nx2 > 1) {
        sm = std::abs(std::sin(pmb->pcoord->x2f(j  )));
        sp = std::abs(std::sin(pmb->pcoord->x2f(j+1)));
        cm = std::cos(pmb->pcoord->x2f(j  ));
        cp = std::cos(pmb->pcoord->x2f(j+1));
      } else {
        // js->jl???
        sm = std::abs(std::sin(pmb->pcoord->x2f(pmb->js  )));
        sp = std::abs(std::sin(pmb->pcoord->x2f(pmb->js+1)));
        cm = std::cos(pmb->pcoord->x2f(pmb->js  ));
        cp = std::cos(pmb->pcoord->x2f(pmb->js+1));
      }
      x2 = pmb->pcoord->x2v(j);
      sinx2 = std::sin(x2);
      src1_j = (sp - sm)/std::abs(cm - cp);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        x1 = pmb->pcoord->x1v(i);
        rm = pmb->pcoord->x1f(i  );
        rp = pmb->pcoord->x1f(i+1);
        src1_i = 1.5*(rp*rp - rm*rm)/(rp*rp*rp - rm*rm*rm);
        Real rho = prim(IDN,k,j,i);
        // Calculate Fluxes
        if (pmb->block_size.nx3 > 1) {
          flux_kp = pmb->phydro->flux[X3DIR](IDN,k+1,j,i);
          flux_km = pmb->phydro->flux[X3DIR](IDN,k,j,i);
        } else {
          flux_kp = rho * prim(IVZ,k,j,i);
          flux_km = rho * prim(IVZ,k,j,i);
        }
        if (pmb->block_size.nx2 > 1) {
          flux_jp = pmb->phydro->flux[X2DIR](IDN,k,j+1,i);
          flux_jm = pmb->phydro->flux[X2DIR](IDN,k,j,i);
        } else {
          flux_jp = rho * prim(IVY,k,j,i);
          flux_jm = rho * prim(IVY,k,j,i);
        }
        flux_ip = pmb->phydro->flux[X1DIR](IDN,k,j,i+1);
        flux_im = pmb->phydro->flux[X1DIR](IDN,k,j,i);
        Real dV = (1./3.)*(rp*rp*rp - rm*rm*rm)*(cm - cp);
        // Calculate Source Terms
        Real cor_1 = omega_t*x1*sinx2*(flux_kp + flux_km)*src1_i;
        Real cen_1 = rho*SQR(omega_t*x1*sinx2)*src1_i;
        Real cor_2 = cor_1*src1_j;
        Real cen_2 = cen_1*src1_j;
        Real cor_3 = -omega_t*sinx2*(cm-cp)*
                    ((std::pow(rp,4) - SQR(rp*x1))*flux_ip +
                     (SQR(rm*x1) - std::pow(rm,4))*flux_im)/(x1*dV);
        Real cen_3 = -omega_t*x1*(rp*rp - rm*rm)*
                     (flux_jp*(sp*sp*sp - sp*SQR(sinx2)) +
                      flux_jm*(sm*SQR(sinx2) - sm*sm*sm))/(2*sinx2*dV);
	pmb->user_out_var(5,k,j,i)=cor_1;
	pmb->user_out_var(6,k,j,i)=cor_2;
	pmb->user_out_var(7,k,j,i)=cor_3;
	pmb->user_out_var(8,k,j,i)=cen_1;
	pmb->user_out_var(9,k,j,i)=cen_2;
	pmb->user_out_var(10,k,j,i)=cen_3;
	cor_[0] = - (cor_2 * cons(IM2,k,j,i) + cor_3 * cons(IM3,k,j,i)) / cons(IM1,k,j,i);
        cor_[1] = - (cor_1 * cons(IM1,k,j,i) + cor_3 * cons(IM3,k,j,i)) / cons(IM2,k,j,i);
        cor_[2] = - (cor_1 * cons(IM1,k,j,i) + cor_2 * cons(IM2,k,j,i)) / cons(IM3,k,j,i);
        if (source_choice == 1) {
          cor_1 = cor_[0];
	}
	else if (source_choice == 3) {
          cor_3 = cor_[2];
	}
	//std::cout<<"cor.mom:   "<<cor_1 * cons(IM1,k,j,i) + cor_2 * cons(IM2,k,j,i)+ cor_3 * cons(IM3,k,j,i)<<"    "<<cor_1<<"    "<<cor_2<<"    "<<cor_3<<"\n";
        //
        cons(IM1,k,j,i) += dt * (cor_1 + cen_1);
        cons(IM2,k,j,i) += dt * (cor_2 + cen_2);
        cons(IM3,k,j,i) += dt * (cor_3 + cen_3);
        // Energy
        if (NON_BAROTROPIC_EOS) {
	  Real e1 = dt*cen_1 * 0.5*(flux_ip + flux_im) / rho;
	  Real e2 = dt*cen_2 * 0.5*(flux_jp + flux_jm) / rho;
          cons(IEN,k,j,i) += e1;
	  cons(IEN,k,j,i) += e2;
	  pmb->user_out_var(11,k,j,i) = e1;
	  pmb->user_out_var(12,k,j,i) = e2;
        }
	//Euler terms
        if (time>=t_L_0) {
          Real vol = pmb->pcoord->GetCellVolume(k,j,i);     
          Real dx3 = pmb->pcoord->dx3f(k);
          Real thr = pmb->pcoord->x2f(j+1);
          Real thl = pmb->pcoord->x2f(j);
          Real eul_src = (-1.0/8.0)*Omega_dot*rho*(dx3/vol)*(std::pow(pmb->pcoord->x1f(i+1),4) - std::pow(pmb->pcoord->x1f(i),4))*(thr-thl-std::sin(thr)*std::cos(thr)+std::sin(thl)*std::cos(thl));
          cons(IM3,k,j,i) += dt * (eul_src);
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i)+= dt*eul_src*0.5*(flux_kp+flux_km)/rho;    
          } 
        }

	pmb->user_out_var(13,k,j,i) = omega_t;
      }
    }
  }
  return;
}

Real grid_space_X2(Real x, RegionSize rs){ //x is the logical location; x=i/nx1, real in [0., 1.]
    Real w=-1.3*x*x*x+1.95*x*x+0.35*x;
    return w*rs.x2max+(1.0-w)*rs.x2min;
}
void UserSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                    const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                    AthenaArray<Real> &cons_scalar) {
  if (rot_flag){
    if (rot_option==2){
      RotFrame2(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
    }
    if (rot_option==1){
      RotFrame(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
    }
      
  }
  if (GENERAL_EOS) {
    //std::cout<<"YES ENTER GENERAL_EOS, include heat_cool term \n";
    heat_cool(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }
}

void HSEInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                FaceField &b, Real time, Real dt, int is, int ie, int js,
                int je, int ks, int ke, int ngh) {
  
  AthenaArray<Real> Ye_T_NR;
  Ye_T_NR.NewAthenaArray(2);  
  Real omega_t = Omega_0;
   

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
        }
      }
    }
  }


  Real kbol = 1.380649e-16;
  Real mn = 1.6726e-24;
  Real r0 = r_0;
  for (int k=ks; k<=ke; ++k) {
    Real phi = pco->x3v(k);
    for (int j=js; j<=je; ++j) {
      //std::cout<<"In inner BC theta ghost values  "<<pco->x2v(js)<<"    "<<pco->x2v(js-1)<<"    "<<pco->x2v(js-2)<<"    "<<pco->x3v(ks)<<"\n";
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        Real r = pco->x1v(is-i);
        Real vphi = omega_t * r * (std::pow(r_0/r,2));
        Real l = r*vphi;


 
        if(vth_choice ==0) {
           if(time==0.0 && Globals::my_rank == 0 && pmb->lid==0 && (j==js) && (k==ks) && (i==1)) {std::cout<<"Setting vth to 0 at inner boundary \n";}
           prim(IVY,k,j,is-i) = 0.0;
        }
 
        if(vth_choice ==1) {
           if(time==0.0 && Globals::my_rank == 0 && pmb->lid==0 && (j==js) && (k==ks) && (i==1)) {std::cout<<"Setting vth to float at inner boundary \n";}
           prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        }
	  

        if (NSCALARS > 0 && ye_index >= 0 && pmb->lid==0) {
        if(time==0.0 && Globals::my_rank == 0 && pmb->lid==0 && (j==js) && (k==ks) && (i==1) && (pmb->loc.lx1==0)) {std::cout<<"Double NR BC at one cell and constant temp BC for all other cells \n";}

           if (time>0.0 and temp_index>=0) {
              GuessTemp_Bc = 8.6173e-11*pmb->pscalars->r(temp_index,k,j,is-i);
              GuessYe_Bc = pmb->pscalars->r(ye_index,k,j,is-i);
              if(GuessTemp_Bc == 0.0) {GuessTemp_Bc = 4.0;}
              if(GuessYe_Bc == 0.0) {GuessYe_Bc = 0.05;}
              //std::cout<<"In inner BC   "<<GuessTemp_Bc<<"    "<<GuessYe_Bc<<"\n";
           }

           if ((Globals::my_rank == 0) && (pmb->lid == 0) && (j==js) && (k==ks) && (i==1) && (pmb->loc.lx1==0)) { //double NR once for each node if only lid = 0 condition is enforced as lid = 0 is satisfied for one meshblock in each node. Additional condition of (Globals::my_rank == 0) is required to do 2D NR only once. 2D NR should be done only once to ensure constant T in both radial and theta direction in the ghost zones.  So changed MPI communication. Correction after plasmoids paper 4 

              //std::cout<<"In inner BC 2D NR    "<<js<<"     "<<je<<"      "<<pco->x1v(is-i)<<"     "<<pco->x2v(j)<<"     "<<pmb->gid<<"\n";
              BC_T0(GuessYe_Bc,GuessTemp_Bc, time, r, theta, Ye_T_NR);
              fixed_tempBC = Ye_T_NR(0);  

           }

           if (n_bc_ranks) {
              if (pmb->lid == 0 && pmb->block_size.x1min == r_0 && (j==js) && (k==ks) && (i==1)) {
                 MPI_Bcast(&fixed_tempBC, 1, MPI_ATHENA_REAL, 0, bc_comm);
              }
           }
        }
        if (NSCALARS > 0 && ye_index >= 0) { 
           if ((Globals::my_rank != 0) || (Globals::my_rank == 0 && pmb->lid!=0) || (Globals::my_rank == 0 && i!=1) || (Globals::my_rank == 0 && i==1 && j!=js)) {

              //std::cout<<"ENTER single NR else    "<<Globals::my_rank<<"    "<<pmb->lid<<"\n";
              BC_T0_gh2(GuessYe_Bc,fixed_tempBC, time, r, theta, Ye_T_NR);
              //std::cout<<"EXIT single NR else    "<<Globals::my_rank<<"    "<<pmb->lid<<"        "<<fixed_tempBC<<"\n";
           }


           T_0=Ye_T_NR(0);
           T_0= T_0/8.6173e-11;

           
           pmb->pscalars->r(ye_index,k,j,is-i) = Ye_T_NR(1);
           dpdd_0 = (kbol*T_0)/mn;
         } 
         else {
           pmb->pscalars->r(0,k,j,is-i) = Ye;
           T_0=zeroQ_temp(0.0, pmb->pscalars->r(0,k,j,is-i),time);
           dpdd_0 = (kbol*T_0)/mn;
           
          }
        
        if (NSCALARS > 0 && temp_index >= 0) {
           pmb->pscalars->r(temp_index,k,j,is-i) = T_0;

         }

	if (GENERAL_EOS) {
	  prim(IDN,k,j,is-i) = rho_0*std::exp((r0 - r) * mu / (dpdd_0 * r * r0)) * std::exp((0.5*omega_t*omega_t/dpdd_0)*(r*r*std::sin(theta)*std::sin(theta)-r0sq));

	  
	}


        if(vr_choice == 0) {
           if(time==0.0 && Globals::my_rank == 0 && pmb->lid==0 && (j==js) && (k==ks) && (i==1)) {std::cout<<"Setting vr to 0 at inner boundary \n";}
           prim(IVX,k,j,is-i) = 0.0;
        }

        if(vr_choice == 1) {
           if(time==0.0 && Globals::my_rank == 0 && pmb->lid==0 && (j==js) && (k==ks) && (i==1)) {std::cout<<"Setting vr by constant Mdot at inner boundary \n";}
           prim(IVX,k,j,is-i) = (prim(IDN,k,j,is)/prim(IDN,k,j,is-i))*prim(IVX,k,j,is) * SQR(pco->x1v(is) / pco->x1v(is-i));
        }    
        if(vr_choice == 2) {
           if(time==0.0 && Globals::my_rank == 0 && pmb->lid==0 && (j==js) && (k==ks) && (i==1)) {std::cout<<"Setting vr to float at inner boundary \n";}
           prim(IVX,k,j,is-i) = prim(IVX,k,j,is);
        }

    if (rot_flag) {
        prim(IVZ,k,j,is-i) = 0.0;
     }
    else {
        prim(IVZ,k,j,is-i) = omega_t*r*std::sin(theta);
    }
        

         if (NON_BAROTROPIC_EOS) {

         Real r_cell[NSCALARS];
         for (int n=0; n<NSCALARS; ++n) {
            r_cell[n] = pmb->pscalars->r(n,k,j,is-i);

         }

            prim(IPR,k,j,is-i) = pmb->peos->PresFromRhoT(prim(IDN,k,j,is-i), T_0,r_cell);
        }
      }

    }
  }
  //std::cout<<"EXITING inner BC   "<<Globals::my_rank<<pmb->lid<<"\n";
}


// Outflow Boundary Condition
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {
  Real omega_t = Omega_0;

  Real re = pco->x1v(ie);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real sinx2 = std::sin(pco->x2v(j));
      for (int i=1; i<=ngh; ++i) {
        Real rgh = pco->x1v(ie+i);
        Real ratio = re / rgh;
        prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie) * SQR(ratio);
        prim(IVX,k,j,ie+i) = std::max(prim(IVX,k,j,ie), 0.0);
        prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie) * ratio;
        prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie) * ratio;
        if (rot_flag) {
	        prim(IVZ,k,j,ie+i) += omega_t*sinx2*re*(ratio - 1.0/ratio);
        }
        if (ye_index>=0){
          pmb->pscalars->r(ye_index,k,j,ie+i)=pmb->pscalars->r(ye_index,k,j,ie);
        }
        else {
          pmb->pscalars->r(0,k,j,ie+i)=Ye;
        }
        if (temp_index>=0){
          pmb->pscalars->r(temp_index,k,j,ie+i)=pmb->pscalars->r(temp_index,k,j,ie);
        }
        if (NON_BAROTROPIC_EOS) {
         Real r_cell[NSCALARS];
         for (int n=0; n<NSCALARS; ++n) {
            r_cell[n] = pmb->pscalars->r(n,k,j,ie+i);
         }
	 Real tempc = pmb->peos->TFromRhoP(prim(IDN,k,j,ie),prim(IPR,k,j,ie),r_cell);
	 prim(IPR,k,j,ie+i) = pmb->peos->PresFromRhoT(prim(IDN,k,j,ie+i),tempc,r_cell);
         
	}
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
        }
      }
    }
  }
}


void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  Real omega_t = Omega_0;


  AthenaArray<Real> le;
  le.NewAthenaArray(12);
  lum_func(pmy_mesh->time,le);
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
	   if (GENERAL_EOS){
         Real r_cell[NSCALARS];
         for (int n=0; n<NSCALARS; ++n) {
            r_cell[n] = pscalars->r(n,k,j,i);
         }
      if(rot_flag){
          ;
      }
      else {
         user_out_var(9,k,j,i) = omega_t;
      }


	  Real temp=peos->TFromRhoP(phydro->w(IDN,k,j,i),phydro->w(IPR,k,j,i),r_cell);
	 

	  user_out_var(0,k,j,i) = temp;
	  Real r = pcoord->x1v(i);

	  Real x=std::pow((1.0-(r_0*r_0)/(r*r)),0.5);
          if (NSCALARS > 0 && ye_index >= 0) {
             Ye = pscalars->r(ye_index,k,j,i);
          }
      Real Eta = QWEta(phydro->w(IDN,k,j,i),temp,Ye);
      temp = temp * 8.6173e-11;  //T in MeV

	  Real qdot = qdotQW(temp, phydro->w(IDN,k,j,i), x, Ye, pmy_mesh->time,Eta,le,0) * 1.6022e-6; // in ergs/s/g


	  user_out_var(1,k,j,i) = qdot;
	  user_out_var(2,k,j,i) = std::sqrt(peos->AsqFromRhoP(phydro->w(IDN,k,j,i),phydro->w(IPR,k,j,i),r_cell));


  Real lnu_y = le(0);
  Real lnubar_y = le(1);
  Real lmu_y = le(2);
  Real eavg1_mu = le(5);
  Real eavg1_nu = le(6) ;
  Real eavg1_nubar = le(8);

	  user_out_var(3,k,j,i) = lnu_y;
	  user_out_var(4,k,j,i) = lnubar_y;
	  user_out_var(5,k,j,i) = lmu_y;
	  user_out_var(6,k,j,i) = eavg1_nu;
	  user_out_var(7,k,j,i) = eavg1_nubar;
	  user_out_var(8,k,j,i) = eavg1_mu;
	}
      }
    }
  }

  if (COMP_DIV_B & MAGNETIC_FIELDS_ENABLED) {
    AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
    int ncells1 = block_size.nx1+2*NGHOST;
    x1area.NewAthenaArray(ncells1+1);
    x2area.NewAthenaArray(ncells1);
    x2area_p1.NewAthenaArray(ncells1);
    x3area.NewAthenaArray(ncells1);
    x3area_p1.NewAthenaArray(ncells1);
    vol.NewAthenaArray(ncells1);
    dflx.NewAthenaArray(ncells1);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        // calculate x1-flux divergence
        pcoord->Face1Area(k,j,is,ie+1,x1area);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          dflx(i) = (x1area(i+1)*pfield->b.x1f(k,j,i+1) - x1area(i)*pfield->b.x1f(k,j,i));
        }

        // calculate x2-flux divergence
        if (block_size.nx2 > 1) {
          pcoord->Face2Area(k,j  ,is,ie,x2area   );
          pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(i) += (x2area_p1(i)*pfield->b.x2f(k,j+1,i)
                       -x2area(i)*pfield->b.x2f(k,j,i));
          }
        }

        // calculate x3-flux divergence
        if (block_size.nx3 > 1) {
          pcoord->Face3Area(k  ,j,is,ie,x3area   );
          pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(i) += (x3area_p1(i)*pfield->b.x3f(k+1,j,i)
                       -x3area(i)*pfield->b.x3f(k,j,i));
          }
        }

        // update conserved variables
        pcoord->CellVolume(k,j,is,ie,vol);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          user_out_var(10,k,j,i) = dflx(i)/vol(i);
        }
      }
    }
    x1area.DeleteAthenaArray();
    x2area.DeleteAthenaArray();
    x2area_p1.DeleteAthenaArray();
    x3area.DeleteAthenaArray();
    x3area_p1.DeleteAthenaArray();
    vol.DeleteAthenaArray();
    dflx.DeleteAthenaArray();
  }
}


void MeshBlock::UserWorkInLoop() {
  AthenaArray<Real>& Omega_pl{pmy_mesh->ruser_mesh_data[5]};
  AthenaArray<Real>& rot_pl{pmy_mesh->ruser_mesh_data[6]};
  if(rot_flag){
    rot_pl(0)=1.0;
  }
  else{
    rot_pl(0)=0.0;
  }

  if (pmy_mesh->time >= t_L_0) {
    int& index_compute{iuser_meshblock_data[0](0)};
    AthenaArray<Real>& my_data = pmy_mesh->ruser_mesh_data[4];

    Real &local_sum{my_data(0)};
    // Calculate cell-centered magnetic field
    AthenaArray<Real> bb;
    bb.NewAthenaArray(3, ke+1, je+1, ie+NGHOST+1);
    pfield->CalculateCellCenteredField(pfield->b, bb, pcoord, is-NGHOST, ie+NGHOST, js,
                                       je, ks, ke);

   
    if (lid == 0) {  // if zeroth block
      local_sum = 0.0;
    }
    if (index_compute >= 0) {
      int i = index_compute;
      Real r = pcoord->x1v(i);
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          Real theta = pcoord->x2v(j);
          Real sth = std::sin(theta);
          Real dcos = std::abs(std::cos(pcoord->x2f(j+1)) - std::cos(pcoord->x2f(j)));
          Real dphi = pcoord->x3f(k+1) - pcoord->x3f(k);

          Real& bcc1 = bb(IB1,k,j,i);
          Real& bcc2 = bb(IB2,k,j,i);
          Real& bcc3 = bb(IB3,k,j,i);
          if(rot_flag){
             local_sum += (phydro->w(IDN,k,j,i)*phydro->w(IVX,k,j,i)*(phydro->w(IVZ,k,j,i)+r*sth*Omega_0)-bcc1*bcc3)*r*r*r*sth*dcos*dphi;
          }
          else {
             local_sum += (phydro->w(IDN,k,j,i)*phydro->w(IVX,k,j,i)*(phydro->w(IVZ,k,j,i))-bcc1*bcc3)*r*r*r*sth*dcos*dphi;
          }
        }
      }
    }
    if (lid == pmy_mesh->nblocal-1) {  // if last block
      Real global_integral;
      MPI_Allreduce(&local_sum, &global_integral, 1,
                  MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      // Set Omega and Omega-dot

      Real MOI = (2.0/5.0)*Mpns*1.9891e33*r_0*r_0;
      Omega_dot = -1.0*global_integral/MOI;
      Omega_0 = Omega_0 + Omega_dot*pmy_mesh->dt;

    }

  }
  Omega_pl(0)=Omega_0;
  //std::cout<<"Omega_pl    "<<Omega_pl(0)<<"\n";
  //std::cout<<"rot_pl    "<<rot_pl(0)<<"\n";
}



















