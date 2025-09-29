#define MYDEBUG1
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file eos_qw.cpp
//! \brief implements functions in class EquationOfState for simple hydrogen EOS
//======================================================================================

// C headers

// C++ headers
#include <algorithm>
#include <cmath>   // sqrt()
#include <fstream>
#include <iostream> // ifstream
#include <limits>   // std::numeric_limits<float>::epsilon()
#include <sstream>
#include <stdexcept> // std::invalid_argument
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../parameter_input.hpp"
#include "../eos.hpp"

// Same effect as adding static to everything inside
namespace {
  // const Real float_eps = std::numeric_limits<Real>::epsilon();
  // const Real float_1pe = 1.0 + float_eps;
  const Real float_max = std::numeric_limits<Real>::max();
  Real prec = 1e-8;
  Real T_floor, T_ceil, LastTemp;
  Real fixed_ye = -1.0;
  int i_ye = -1;
  int i_temp = -1;
  bool use_T_floor;
  int nmax;
  AthenaArray<Real> EosData;

  const Real third = 1.0 / 3.0;
  const Real c = 2.99792458e10;
  const Real k = 1.380649e-16;
  const Real mn = 1.6726e-24;
  const Real hbar = 6.62607015e-27/(2.0*PI);
  const Real c3 = std::pow(k/(hbar*c),3);
  const Real con3 = (11.0*PI*PI/180.0)*c3*k;
  const Real Tmin=1e3;
  const Real Tmax=1e15;
  const Real kmn=k/mn;
  const Real thel= 30.0/11.0;
  const Real fiel= 15.0/11.0;
  const Real eta_den_const = std::pow(6.0,2.0*third);
} // namespace

// EOS data indicies
namespace EOS {
  enum EosIndex {iE=0, idEdT=1, iP=2, idPdT=3, iAsq=4, iT=5, N=6};
}


Real PorE(Real rho, Real T, Real Ye, int index) {
  Real vol = mn/rho;
  Real T4 = T*T*T*T;
  Real a = c3*std::pow(T,3.0)*vol*(PI/3.0);
  Real a2 = SQR(a);
  Real a4 = SQR(a2);
  Real a6 = a2 * a4;
  Real y2 = SQR(Ye);
  Real b = std::sqrt(4.0*a6+27.0*a4*y2);
  Real term = std::pow(9.0*a2*Ye+std::pow(3.0,0.5)*b, third);
  Real eta = (std::pow(2.0,third)/eta_den_const)*term/a - (2.0*std::pow(3.0,third)/eta_den_const)*a/term; // actually eta/pi

  Real eta2 = SQR(eta);

  Real eta4 = SQR(eta2);

  Real con= con3*(1.0+thel*eta2 + fiel*eta4);
  Real p0 = rho*kmn*T;
  if (index == EOS::iP) {
    return con*T4 + p0;
  } else if (index == EOS::iE) {
    return 3.0*con*T4 + 1.5*p0;
  }
  std::cout<<"Incorrect index in PorE \n";
  return 0;
}




void QWData(Real rho, Real T, Real Ye, AthenaArray<Real> &OutData) {
  //if (use_T_floor){
    //std::cout<<"ENTER IF IN QWDATA \n";
    //T = std::max(T, T_floor);
  //}
  Real vol = mn/rho;
  Real T3 = T*T*T;
  Real T4 = T*T3;
  Real a = c3*std::pow(T,3)*vol*(PI/3.0);
  Real a2 = SQR(a);
  Real a4 = SQR(a2);
  Real a6 = a2 * a4;
  Real y2 = SQR(Ye);
  Real b = std::sqrt(4.0*a6+27.0*a4*y2);
  Real term = std::pow(9.0*a2*Ye+std::pow(3.0,0.5)*b, third);
  Real eta = (std::pow(2.0,third)/eta_den_const)*term/a - (2.0*std::pow(3.0,third)/eta_den_const)*a/term; // actually eta/pi
  
  Real eta2 = SQR(eta);
  Real eta3 = eta*eta2;
  Real eta4 = SQR(eta2);

  Real da_dt = 3.0*a/T;
  Real dterm_dt= third*std::pow(term, -2)*(18.0*a*Ye*da_dt
                + std::pow(3.0,0.5)*0.5*a*(24.0*a4*da_dt+27.0*4.0*a2*y2*da_dt)/b);
  Real der_eta = (std::pow(2.0,third)/eta_den_const)*(dterm_dt/a - term*da_dt/(a*a))
                 - (2.0*std::pow(3.0,third)/eta_den_const)*(da_dt/term -a*dterm_dt/(term*term)); //actually derivative of eta/PI

  Real con= con3*(1.0+thel*eta2 + fiel*eta4);
  Real p0 = rho*kmn*T;
  Real P = con*T4 + p0;
  Real e = 3.0*con*T4 + 1.5*p0;

  Real der_P= 4.0*T3*con + con3*T4*(thel*2.0*eta*der_eta
                                    + fiel*4.0*eta3*der_eta);
  Real der_e= 3.0*der_P+ 1.5*rho*kmn;
  der_P += rho*kmn;

  Real deta_drho = (eta/rho)*(1.0+eta2)/(1.0+3.0*eta2); //actually rho derivative of eta/PI
  Real ctsq=kmn*T + con3*T4*4.0*fiel*eta*deta_drho*(1.0+eta2);
  Real D=(T/rho)*der_P;
  Real cv= der_e / rho;
  Real asq = ctsq+D*D/(cv*T);

  OutData(0) = e;
  OutData(1) = der_e;
  OutData(2) = P;
  OutData(3) = der_P;
  OutData(4) = asq;
  OutData(5) = T;

#ifdef MYDEBUG1
  bool flag = false;
  for (int i=0; i<6; ++i) {
    if (!std::isfinite(OutData(i)))
      flag = true;
  }
  if (flag) {
    printf("ERR: Non-finite output detected in EOS output.\n");
    printf("e_int, de/dT, P_gas, dP/dT,   a^2,  Temp,  dens\n");
    printf("EOS data: ");
    for (int i=0; i<6; ++i)
      printf("%.4e, ", OutData(i));
    printf("%.4e\n", rho);
    std::stringstream msg;
    msg << "### FATAL ERROR in QWData" << std::endl
        << "Non-finite value in EOS data." << std::endl;
    ATHENA_ERROR(msg);
  }
#endif // MYDEBUG1
}




// index = 0 for internal energy; index = 2 for pressure; var = int energy or pressure
void TempInvert(Real rho, Real GuessTemp, Real Ye, Real var, const int index,
                AthenaArray<Real> &OutData) {
  if ((std::isnan(var)) || (std::isnan(rho)) || (std::isnan(GuessTemp))) {
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("ERR (%s): %.4e, rho: %.4e,  Temp: %.4e \n", varnames[index],var,rho,
             GuessTemp);

      std::stringstream msg;
      msg <<"Nan in root find (var) \n";
      ATHENA_ERROR(msg);
    }
  if ((var<=0.0) || (rho<=0.0)) {
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("ERR (%s): %.4e, rho: %.4e,  Temp: %.4e \n", varnames[index],var,rho,
             GuessTemp);

      std::stringstream msg;
      msg <<"Negative var\n";
      ATHENA_ERROR(msg);
  }


  Real BrakT[] = {T_floor,T_ceil};
  Real BrakVal[] = {0, 0};
  Real InvVar = 1.0 / var;
  Real error = float_max;
  int nlim = nmax;
  Real LastTemp = BrakT[0];
  QWData(rho, LastTemp,Ye, OutData);
  BrakVal[0] = OutData(index) * InvVar - 1.0;
  Real LastErr = BrakVal[0];
  Real delta;
#ifdef MYDEBUG1
  //printf("%d: %.16e, %.16e\n", index, var, rho);
  int mode = 0;
#endif
  while (std::abs(error) > prec) {
    if(std::isnan(GuessTemp)) {

      std::stringstream msg;
      msg <<"Nan in GuessTemp in root find loop \n";
      ATHENA_ERROR(msg);
    }
    if (BrakVal[0] > 0) {
      QWData(rho, BrakT[0],Ye, OutData);
      Real low = OutData(index);
      // If we are below Tmin
      if (var < low) {
        Real rtemp=1e3;
        QWData(rho, rtemp,Ye, OutData);
        
        return;
      }
      QWData(rho, BrakT[1],Ye, OutData);
      Real high =  OutData(index);
      std::stringstream msg;
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("0 ERR (%s): %.4e !<= %.4e !<= %.4e\n", varnames[index], low, var,
             high);
      printf("at rho = %.4e, T_bounds = %.4e, %.4e,\n", rho, BrakT[0], BrakT[1]);

      msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
          << std::endl << "Root not bracketed" << std::endl;
      ATHENA_ERROR(msg);
    }

    // if we are outside brackets use bisection method for a step
    if ((GuessTemp <= BrakT[0]) || (GuessTemp >= BrakT[1])) {
      GuessTemp = std::sqrt(BrakT[0] * BrakT[1]);

#ifdef MYDEBUG1
      mode = 1;
#endif
    }
    QWData(rho, GuessTemp,Ye, OutData);
    error = OutData(index) * InvVar - 1.0;
#ifdef MYDEBUG1
    if(GuessTemp<0.0 ||(std::isnan(GuessTemp))) {
      printf("%04d [%.4g, %.4g, %.4g]; %.4g| %d\n", 10000 - nlim, BrakT[0], GuessTemp,
             BrakT[1], error, mode);
    }
#endif
    // update bracketing values
    if (error < 0) {
      BrakT[0] = GuessTemp;
      BrakVal[0] = error;
    } else {
      BrakT[1] = GuessTemp;
      BrakVal[1] = error;
    }
    if (BrakT[1] <= BrakT[0]) {
      QWData(rho, BrakT[0],Ye, OutData);
      Real low = OutData(index);
      // If we've specified Tfloor and we are below Tmin just use Tmin and return
      if (var < low) {
        Real rtemp=1e3;
        QWData(rho, rtemp,Ye, OutData);

        return;
      }
      QWData(rho, BrakT[1],Ye, OutData);
      Real high = OutData(index);
      std::stringstream msg;
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("1 ERR (%s): %.4e !<= %.4e !<= %.4e\n", varnames[index], low, var,
             high);
      printf("at rho = %.4e, T_bounds = %.4e, %.4e,\n", rho, BrakT[0], BrakT[1]);

      msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
          << std::endl << "Root not bracketed" << std::endl;
      ATHENA_ERROR(msg);
    }

    if (std::abs(error) > 1e-2) {
      // secant method step
      delta = error * (GuessTemp - LastTemp) / (error - LastErr);
      LastTemp = GuessTemp;
      GuessTemp -= delta;
      //if(std::isnan(GuessTemp)){
      //  std::cout<<"Nan in secant   "<<rho<<"     "<<delta<<"    "<<error-LastErr
      //           <<"    "<<GuessTemp<<"    "<<LastTemp<<"    "<<error<<"\n";
      //}
      LastErr = error;
#ifdef MYDEBUG1
      mode = 2;
#endif
    } else {
      // Newton–Raphson step
      delta = var * error / OutData(index + 1);
      GuessTemp -= delta;
      //if(std::isnan(GuessTemp)){
      //  std::cout<<"Nan in NR  "<<OutData(index+1)<<"    "<<delta<<"\n";
      //}
#ifdef MYDEBUG1
      mode = 3;
#endif
    }

    if (nlim-- < 0) {
      
      std::stringstream msg;
      msg <<"Nlim reached in TempInvert  \n";
      ATHENA_ERROR(msg);
      return;
    }
  }


  //if(nmax-nlim>100){std::cout<<"Iter     "<<nmax-nlim<<"\n";}
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//  \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas, Real* s) {
  Real ye = fixed_ye;

  if (NSCALARS > 0 && i_ye >= 0) {

      ye = s[i_ye]/rho;
   
   }
  if (NSCALARS > 0 && i_temp >= 0) {
    LastTemp = s[i_temp]/rho;
  }
  TempInvert(rho, LastTemp, ye, egas, EOS::iE, EosData);
  LastTemp = EosData(EOS::iT);
  if (NSCALARS > 0 && i_temp >= 0) {
    s[i_temp] = LastTemp * rho;
  }
  return EosData(EOS::iP);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//  \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres, Real* r) {
  Real ye = fixed_ye;
  if (NSCALARS > 0 && i_ye >= 0) {

      ye = r[i_ye];
  }   
  
  if (NSCALARS > 0 && i_temp >= 0) {
    LastTemp = r[i_temp];
  }
  //std::cout<<"In EOS egasfromrhoP   "<<ye<<"\n";
  TempInvert(rho, LastTemp, ye, pres, EOS::iP, EosData);
  LastTemp = EosData(EOS::iT);
  if (NSCALARS > 0 && i_temp >= 0) {
    r[i_temp] = LastTemp;
  }
  return EosData(EOS::iE);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//  \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres, const Real* r) {
  Real ye = fixed_ye;
  if (NSCALARS > 0 && i_ye >= 0) {

      ye = r[i_ye];
  
  }
  if (NSCALARS > 0 && i_temp >= 0) {
    LastTemp = r[i_temp];
  }
  //std::cout<<"In EOS Asq   "<<ye<<"\n";
  TempInvert(rho, LastTemp, ye, pres, EOS::iP, EosData);
  LastTemp = EosData(EOS::iT);

  return EosData(EOS::iAsq);
}


Real EquationOfState::TFromRhoP(Real rho, Real pres, Real* r) {
  Real ye = fixed_ye;
  if (NSCALARS > 0 && i_ye >= 0) {
    ye = r[i_ye];
  }
  if (NSCALARS > 0 && i_temp >= 0) {
    LastTemp = r[i_temp];
  }
  TempInvert(rho, LastTemp, ye, pres, EOS::iP, EosData);
  LastTemp = EosData(EOS::iT);
  if (NSCALARS > 0 && i_temp >= 0) {
    r[i_temp] = LastTemp;
  }
  return LastTemp;
}

Real EquationOfState::TFromRhoEgas(Real rho, Real egas, Real* s) {
  Real ye = fixed_ye;
  if (NSCALARS > 0 && i_ye >= 0) {
    ye = s[i_ye] / rho;
  }

  if (NSCALARS > 0 && i_temp >= 0) {
    LastTemp = s[i_temp]/rho;
  }

  TempInvert(rho, LastTemp, ye, egas, EOS::iE, EosData);
  LastTemp = EosData(EOS::iT);
  if (NSCALARS > 0 && i_temp >= 0) {
    s[i_temp] = LastTemp*rho;
  }
  return LastTemp;
}

Real EquationOfState::PresFromRhoT(Real rho, Real T, Real* r) {
  Real ye = fixed_ye;
  if (NSCALARS > 0 && i_ye >= 0) {
    ye = r[i_ye];
  }

  QWData(rho, T, ye, EosData);
  return EosData(EOS::iP);
}

Real EquationOfState::GetEgasFloor(Real rho, Real* s) {
  Real ye = fixed_ye;
  if (NSCALARS > 0 && i_ye >= 0) {
    ye = s[i_ye] / rho;
  }
  QWData(rho, T_floor, ye, EosData);
  return EosData(EOS::iE);
  
}

Real EquationOfState::GetPresFloor(Real rho, Real* r) {
  Real ye = fixed_ye;
  if (NSCALARS > 0 && i_ye >= 0) {
    ye = r[i_ye];
  }
  QWData(rho, T_floor, ye, EosData);
  return EosData(EOS::iP);
  
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput* pin) {
  EosData.NewAthenaArray(EOS::N);
  prec = pin->GetOrAddReal("hydro", "InversionPrecision", prec);
  if (pin->DoesParameterExist("hydro", "eos_ye")) {
    std::cout<<"Using fixed Ye in QW EOS \n";
    fixed_ye = pin->GetReal("hydro", "eos_ye");
  }
  if (pin->DoesParameterExist("hydro", "eos_ye_index")) {
    if (Globals::my_rank == 0) {std::cout<<"Using Ye passive scalars in QW EOS \n";}
    i_ye = pin->GetInteger("hydro", "eos_ye_index");
    if (i_ye < 0 || i_ye >= NSCALARS) {
      std::stringstream msg;
      msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
          << "hydro/eos_ye_index must be between 0 and NSCALARS (" << NSCALARS << ")."
          << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  if (fixed_ye < 0 && i_ye < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
        << "either hydro/eos_ye or hydro/eos_ye_index must be specified."
        << std::endl;
    ATHENA_ERROR(msg);
  }

    if (pin->DoesParameterExist("hydro", "helm_temp_index")) {
      i_temp = pin->GetInteger("hydro", "helm_temp_index");
      if (Globals::my_rank == 0) {std::cout<<"Using Temp passive scalars in EOS\n";}
      if (i_temp < 0 || i_temp >= NSCALARS) {
        std::stringstream msg;
        msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
            << "hydro/helm_temp_index must be between 0 and NSCALARS ("
            << NSCALARS << ")."
            << std::endl;
        ATHENA_ERROR(msg);
      }
    }

  T_floor = pin->GetOrAddReal("hydro", "T_floor", Tmin);
  T_ceil = pin->GetOrAddReal("hydro", "T_ceil", Tmax);
  if(T_floor<Tmin){
    T_floor=Tmin;
  }
  if (Globals::my_rank == 0) {std::cout<<"T floor and ceil   "<<T_floor<<"   "<<T_ceil<<"\n";}
  use_T_floor = pin->GetOrAddBoolean("hydro", "use_T_floor", true);
  LastTemp = pin->GetOrAddReal("hydro", "T0", std::sqrt(T_floor * T_ceil));
  nmax = pin->GetOrAddInteger("hydro", "nmax", 10000);
  // if root processor and zeroth block
  if ((Globals::my_rank == 0) && (pmy_block_->lid == 0)) {
    //std::cout<<Ye<<" Ye\n";
    std::cout<<"Using updated EOS QW NEW 3 \n";
  }
  #ifdef MYDEBUG
    Real rho, temp;
    std::cout << "Input fluid parameters and retrieve EOS parameters." << '\n'
              << "Non-positive inputs will exit loop." << '\n';
    std::cout << "Input density (mass/volume): ";
    std::cin >> rho;
    std::cout << "Input temperature (K): ";
    std::cin >> temp;

    while (rho > 0 && std::isfinite(rho) && temp > 0 && std::isfinite(temp)) {
      printf("d, t: %.16e, %.16e\n", rho, temp);
      QWData(rho, temp, EosData);
      Real p = EosData(EOS::iP);
      Real e = EosData(EOS::iE);
      Real a2 = EosData(EOS::iAsq);
      Real P_tejas = compare::P_of_rho_T(rho, temp);
      Real e_tejas = compare::e_of_rho_T(rho, temp);
      Real Asq_tejas = compare::asq(rho, temp, P_tejas);
      printf("e(d, T)          , p(d, T)         , ASq(d, T)\n");
      printf("%.16e, %.16e, %.16e\n", e, p, a2);
      printf("e_matt/e_tejas-1, P_matt/P_tejas-1, Asq_matt/Asq_tejas-1\n");
      printf("%.16e, %.16e, %.16e\n", e/e_tejas-1.0, p/P_tejas-1.0, a2/Asq_tejas-1.0);
      Real Te = TFromRhoEgas(rho, e);
      Real Tp = TFromRhoP(rho, p);
      printf("T(d, e)/T-1, T(d, p)/T-1\n");
      printf("%.16e, %.16e\n", Te/temp - 1.0, Tp/temp - 1.0);
      std::cout << "Input density (mass/volume): ";
      std::cin >> rho;
      std::cout << "Input temperature (K): ";
      std::cin >> temp;
    }
    std::cout << std::endl;
  #endif
  return;
}




