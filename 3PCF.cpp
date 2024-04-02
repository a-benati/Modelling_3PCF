/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/**
 *  @file Cosmology/Lib/3PCF.cpp
 *
 *  @brief Methods of the class ThreePointCorrelation used to model three-point
 *  statistics
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation, used to model three-point statistics
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "3PCF.h"

using namespace std;

using namespace cbl;
using namespace cosmology;


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::denominator_Q (const double r1, const double r2, const double theta, const vector<double> rr, const vector<double> xi_matter) const
{
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*cos(theta));
  const glob::FuncGrid interp_xi_matter(rr, xi_matter, "Spline");

  const double xi1 = interp_xi_matter(r1);
  const double xi2 = interp_xi_matter(r2);
  const double xi3 = interp_xi_matter(r3);

  return xi1*xi2+xi2*xi3+xi3*xi1;
}


// =====================================================================================


void cbl::cosmology::ThreePointCorrelation::integrals_Q_nonLocal (vector<double> &xi_matter, vector<double> &Phi, const vector<double> rr, const vector<double> kk, const vector<double> Pk_matter, const double prec) const
{
  xi_matter = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_matter, 0, 0, 1, 1);

  const int nk = kk.size();
  glob::FuncGrid interp_Pk = glob::FuncGrid(kk, Pk_matter, "Spline", BinType::_logarithmic_);

  Phi.erase(Phi.begin(), Phi.end());
  Phi.resize(nk, 0);

  for (size_t i=0; i<rr.size(); i++) {

    auto integrand = [&] (const double _k)
    {
      const double kr = _k*rr[i];
      return pow(TopHat_WF(kr), 2)*interp_Pk(_k)*sin(kr)/kr; 
    };
    
    Phi[i] = 1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag(integrand, 0, 100, prec);
  }
  
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::Gamma_3PCF (const double r1, const double r2, const double theta, const vector<double> xi, const vector<double> dPhi) const
{
  return ((xi[0]+3*dPhi[0]/r1)*(xi[1]+3*dPhi[1]/r2))*legendre_polynomial(cos(theta), 2);
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::Q_nonLocal (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  if (rr.size()==0) {
    rr = linear_bin_vector(200, 1., 300.);
    integrals_Q_nonLocal(xi_matter, Phi, rr, kk, Pk_matter, 1.e-3);
  }

  glob::FuncGrid interp_xiDM(rr, xi_matter, "Spline");
  glob::FuncGrid interp_Phi(rr, Phi, "Spline");

  const double mu12 = cos(theta);
      
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*mu12);
  double a12 = theta;

  double a23 = 0.;
  double a13 = 0.; 
  if (r1<=r2) {
    a23 = asin(sin(a12)*r1/r3);
    a13 = par::pi-a12-a23;
  }
  else {
    a13 = asin(sin(a12)*r2/r3);
    a23 = par::pi-a12-a13;
  }

  const double xi1 = interp_xiDM(r1);
  const double xi2 = interp_xiDM(r2);
  const double xi3 = interp_xiDM(r3);
  
  const double dPhi1 = interp_Phi.D1v(r1);
  const double dPhi2 = interp_Phi.D1v(r2);
  const double dPhi3 = interp_Phi.D1v(r3);
  
  double gamma = 0.;
  gamma += Gamma_3PCF(r1, r2, a12, {xi1, xi2}, {dPhi1, dPhi2});
  gamma += Gamma_3PCF(r2, r3, a23, {xi2, xi3}, {dPhi2, dPhi3});
  gamma += Gamma_3PCF(r3, r1, a13, {xi3, xi1}, {dPhi3, dPhi1});
  
  return 2./3*(gamma/denominator_Q(r1, r2, theta, rr, xi_matter)-1);
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::Q_nonLocal (const double r1, const double r2, const std::vector<double> theta, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_matter, Phi;
  vector<double> qNL(ntheta, 0);

  for (int i=0; i<ntheta; i++)
    qNL[i] = Q_nonLocal(r1, r2, theta[i], rr, xi_matter, Phi, kk, Pk_matter);

  return qNL;
}


// =====================================================================================


void cbl::cosmology::ThreePointCorrelation::integrals_zeta_Slepian (std::vector<double> &xi_matter, std::vector<double> &xi_matter_m1, std::vector<double> &xi_matter_p1, std::vector<double> &xi_matter_2, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  vector<double> Pk_matter_m1 = Pk_matter, Pk_matter_p1 = Pk_matter;
  const int nk = kk.size();

  for (int i=0; i< nk; i++) {
    Pk_matter_m1[i] *= pow(kk[i], -1);
    Pk_matter_p1[i] *= kk[i];
  }

  xi_matter = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_matter, 0, 0, 1, 1);
  xi_matter_m1 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_matter_m1, 1, 0, 1, 1);
  xi_matter_p1 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_matter_p1, 1, 0, 1, 1);
  xi_matter_2 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_matter, 2, 0, 1, 1);
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::Q_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &xi_matter_m1, std::vector<double> &xi_matter_p1, std::vector<double> &xi_matter_2, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders, const double prec) const
{
  const double zeta_DM = zeta_DM_Slepian(r1, r2, theta, rr, xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, kk, Pk_matter, norders, prec);

  return zeta_DM/denominator_Q(r1, r2, theta, rr, xi_matter);
}


// =====================================================================================


void cbl::cosmology::ThreePointCorrelation::integrals_zeta_BarrigaGatzanaga (std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int nk = kk.size();

  vector<double> Pk_matter_m4 = Pk_matter;
  for (int i=0; i< nk; i++)
    Pk_matter_m4[i] *= pow(kk[i], -2);
  
  xi_matter = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_matter, 0, 0, 1, 1);
  Phi = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_matter_m4, 0, 0, 1, 1);
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_single_BarrigaGatzanaga (const double r1, const double r2, const double theta, const std::vector<double> xi, const std::vector<double> dxi, const std::vector<double> dPhi) const
{
  const double mu = cos(theta);

  const double xi1 = xi[0];
  const double xi2 = xi[1];

  const double dxi1 = dxi[0];
  const double dxi2 = dxi[1];
 
  const double dPhi1 = dPhi[0];
  const double dPhi2 = dPhi[1];

  const double t1 = (10./7)*xi1*xi2;
  const double t2 = -3*dPhi1*dPhi2/(r1*r2)-xi1*dPhi2/r2-xi2*dPhi1/r1;
  const double t3 = mu*mu*( (xi1+3*dPhi1/r1)*(xi2+3*dPhi2/r2));
  const double t4 = -mu*(dxi1*dPhi2+dxi2*dPhi1);

  return t1+(4./7)*(t2+t3)+t4;
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  if (rr.size()==0) {
    rr = linear_bin_vector(200, 1., 300.);
    integrals_zeta_BarrigaGatzanaga (xi_matter, Phi, rr, kk, Pk_matter);
  }

  glob::FuncGrid interp_xiDM(rr, xi_matter, "Spline");
  glob::FuncGrid interp_Phi(rr, Phi, "Spline");

  const double mu12 = cos(theta);
  
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*mu12);
  const double a12 = theta;
  
  double a23 = 0.;
  double a13 = 0.; 
  if (r1<=r2) {
    a23 = asin(sin(a12)*r1/r3);
    a13 = par::pi-a12-a23;
  }
  else {
    a13 = asin(sin(a12)*r2/r3);
    a23 = par::pi-a12-a13;
  }

  const double xi1 = interp_xiDM(r1);
  const double xi2 = interp_xiDM(r2);
  const double xi3 = interp_xiDM(r3);
  
  const double dxi1 = interp_xiDM.D1v(r1);
  const double dxi2 = interp_xiDM.D1v(r2);
  const double dxi3 = interp_xiDM.D1v(r3);
  
  const double dPhi1 = interp_Phi.D1v(r1);
  const double dPhi2 = interp_Phi.D1v(r2);
  const double dPhi3 = interp_Phi.D1v(r3);

  const double t1 = zeta_single_BarrigaGatzanaga(r1, r2, a12, {xi1, xi2}, {dxi1, dxi2}, {dPhi1, dPhi2});
  const double t2 = zeta_single_BarrigaGatzanaga(r2, r3, a23, {xi2, xi3}, {dxi2, dxi3}, {dPhi2, dPhi3});
  const double t3 = zeta_single_BarrigaGatzanaga(r3, r1, a13, {xi3, xi1}, {dxi3, dxi1}, {dPhi3, dPhi1});

  const double zeta = t1+t2+t3;
  return zeta;
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::Q_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  return zeta_DM_BarrigaGatzanaga(r1, r2, theta, rr, xi_matter, Phi, kk, Pk_matter)/denominator_Q(r1, r2, theta, rr, xi_matter);
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::zeta_DM (const double r1, const double r2, const std::vector<double> theta, const string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_matter;
  vector<double> zDM(ntheta, 0.);

  if (model=="Slepian") {
    vector<double> xi_matter_m1, xi_matter_p1, xi_matter_2;
    for (int i=0; i<ntheta; i++)
      zDM[i] = zeta_DM_Slepian(r1, r2, theta[i], rr, xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, kk, Pk_matter, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga") {
    vector<double> Phi;
    for (int i=0; i<ntheta; i++)
      zDM[i] = zeta_DM_BarrigaGatzanaga(r1, r2, theta[i], rr, xi_matter, Phi, kk, Pk_matter);
  }
  else
    ErrorCBL("the chosen model is not implemented!", "zeta_DM", "3PCF.cpp");

  return zDM;
}

// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::Q_DM (const double r1, const double r2, const std::vector<double> theta, const string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_matter;
  vector<double> qDM(ntheta, 0);

  if (model=="Slepian") {
    vector<double> xi_matter_m1, xi_matter_p1, xi_matter_2;
    for (int i=0; i<ntheta; i++)
      qDM[i] = Q_DM_Slepian(r1, r2, theta[i], rr, xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, kk, Pk_matter, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<ntheta; i++)
      qDM[i] = Q_DM_BarrigaGatzanaga(r1, r2, theta[i], rr, xi_matter, Phi, kk, Pk_matter);
  }
  else
    ErrorCBL("the chosen model is not implemented!", "Q_DM", "3PCF.cpp");

  return qDM;

}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int ntheta = theta.size();
  vector<double> qDM = ThreePointCorrelation::Q_DM(r1, r2, theta, model, kk, Pk_matter);
  vector<double> qH(ntheta, 0);

  for (int i=0; i<ntheta; i++)
    qH[i] = qDM[i]/b1+b2/(b1*b1); 

  return qH;
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double g2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int ntheta = theta.size();
  vector<double> qH = ThreePointCorrelation::Q_halo(r1, r2, theta, b1, b2, model, kk, Pk_matter);
  vector<double> qNL = ThreePointCorrelation::Q_nonLocal(r1, r2, theta, kk, Pk_matter);

  for (int i=0; i<ntheta; i++)
    qH[i] += g2/b1*qNL[i]; 

  return qH;
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::zeta_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int nr = rr.size();
  const double theta = par::pi/3;
  vector<double> _rr, xi_matter;
  vector<double> zDM(nr, 0);

  if (model=="Slepian") {
    vector<double> xi_matter_m1, xi_matter_p1, xi_matter_2;
    for (int i=0; i<nr; i++)
      zDM[i] = zeta_DM_Slepian(rr[i], rr[i], theta, _rr, xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, kk, Pk_matter, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<nr; i++)
      zDM[i] = zeta_DM_BarrigaGatzanaga(rr[i], rr[i], theta, _rr, xi_matter, Phi, kk, Pk_matter);
  }
  else
    ErrorCBL("the chosen model is not implemented!", "zeta_DM_eq", "3PCF.cpp");

  return zDM;
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::Q_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int nr = rr.size();
  const double theta = par::pi/3;
  vector<double> _rr, xi_matter;
  vector<double> Q_DM(nr, 0);

  if (model=="Slepian") {
    vector<double> xi_matter_m1, xi_matter_p1, xi_matter_2;
    for (int i=0; i<nr; i++)
      Q_DM[i] = Q_DM_Slepian(rr[i], rr[i], theta, _rr, xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, kk, Pk_matter, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<nr; i++)
      Q_DM[i] = Q_DM_BarrigaGatzanaga(rr[i], rr[i], theta, _rr, xi_matter, Phi, kk, Pk_matter);
  }
  else
    ErrorCBL("the chosen model is not implemented!", "Q_DM_eq", "3PCF.cpp");

  return Q_DM;
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_multipoles_covariance (const double Volume, const double nObjects, const int l, const int l_prime, const double r1, const double r2, const double r1_prime, const double r2_prime, const double deltaR, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> rr, const std::vector<double> Xi, const double prec)
{
  (void)prec;

  const double kr0 = 1.; //check

  size_t nk = static_cast<int>(kk.size());
  size_t nr = static_cast<int>(rr.size());
  double dlogr = log(rr[1])-log(rr[0]);

  const double density = nObjects/Volume;

  // r binning

  function<double(double, int, double)> func1; 
  function<double(double, int, int, double, double)> func2; 
  function<double(double, double)> func1_noise; 

  if (deltaR<0) {
    func1 = [&] (const double kk, const int l, const double rr) {
      return jl(kk*rr, l); };
    func2 = [&] (const double kk, const int l, const int lp, const double rr, const double rrp) {
      return jl(kk*rr, l)*jl(kk*rrp, lp); };
  } else {
    func1 = [&] (const double kk, const int l, const double rr) {
      return jl_distance_average(kk, l, rr-deltaR*0.5, rr+deltaR*0.5); };
    func2 = [&] (const double kk, const int l, const int lp, const double rr, const double rrp) {
      return jl_distance_average(kk, l, rr-deltaR*0.5, rr+deltaR*0.5)*jl_distance_average(kk, lp, rrp-deltaR*0.5, rrp+deltaR*0.5); };
  }

  func1_noise = [&] (const double rr, const double bin) {
    double rmin = bin-deltaR/2;
    double rmax = bin+deltaR/2;
    if (rr>=rmin && rr<=rmax && deltaR>0)
      return 1./(4*par::pi/3*(pow(rmax, 3)-pow(rmin, 3)));
    else 
      return 0.;
  };

  // integrand 
  
  // Terms of type  I_ll = int [(Pk+1/n) jl] jl
    
  vector<double> pk_r1 = Pk, pk_r2 = Pk;
  vector<double> pk_r1p = Pk, pk_r2p = Pk;

  for (size_t i=0; i<nk; i++) {
    pk_r1[i] *= func1(kk[i], l, r1);
    pk_r2[i] *= func1(kk[i], l, r2);
    pk_r1p[i] *= func1(kk[i], l_prime, r1_prime);
    pk_r2p[i] *= func1(kk[i], l_prime, r2_prime);
  }

  vector<double> I1_r_r1 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1, l, 0, kr0, 1);
  vector<double> I1_r_r2 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2, l, 0, kr0, 1);
  vector<double> I1_r_r1p = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1p, l_prime, 0, kr0, 1);
  vector<double> I1_r_r2p = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2p, l_prime, 0, kr0, 1);

  for (size_t i=0; i<nr; i++) {
    I1_r_r1[i] += func1_noise(rr[i], r1)/density;
    I1_r_r2[i] += func1_noise(rr[i], r2)/density;
    I1_r_r1p[i] += func1_noise(rr[i], r1_prime)/density;
    I1_r_r2p[i] += func1_noise(rr[i], r2_prime)/density;
  }

  // l2
  vector<int> l2;
  for (int ll=abs(l-l_prime); ll<l+l_prime+1; ll++)
    l2.push_back(ll);

  const size_t n_l2 = l2.size();

  // Terms of type  I_l1_l2_l = int [(Pk+1/n) jl1 jl2] jl
  
  vector<double> pk_r1_r1p = Pk, pk_r2_r2p = Pk, pk_r2_r1p = Pk, pk_r1_r2p = Pk;
  for (size_t i=0; i<nk; i++) {
    pk_r1_r1p[i] = (pk_r1_r1p[i]+1./density)*func2(kk[i], l, l_prime, r1, r1_prime);
    pk_r2_r2p[i] = (pk_r2_r2p[i]+1./density)*func2(kk[i], l, l_prime, r2, r2_prime);
    pk_r2_r1p[i] = (pk_r2_r1p[i]+1./density)*func2(kk[i], l, l_prime, r2, r1_prime);
    pk_r1_r2p[i] = (pk_r1_r2p[i]+1./density)*func2(kk[i], l, l_prime, r1, r2_prime);
  }

  vector<vector<double>> I2_r1_r1p(n_l2, vector<double>(nr, 0.));
  vector<vector<double>> I2_r2_r2p(n_l2, vector<double>(nr, 0.));
  vector<vector<double>> I2_r2_r1p(n_l2, vector<double>(nr, 0.));
  vector<vector<double>> I2_r1_r2p(n_l2, vector<double>(nr, 0.));

  for (size_t ll=0; ll<n_l2; ll++) {
    double wig = gsl_sf_coupling_3j(2*l, 2*l_prime, 2*l2[ll], 0, 0, 0);
    if (wig!=0) {
      I2_r1_r1p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1_r1p, l2[ll], 0, kr0, 1);
      I2_r2_r2p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2_r2p, l2[ll], 0, kr0, 1);
      I2_r2_r1p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2_r1p, l2[ll], 0, kr0, 1);
      I2_r1_r2p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1_r2p, l2[ll], 0, kr0, 1);
    }
  }

  // integrand

  double Int = 0;

  for (size_t i=0; i<nr; i++) {
    double sum = 0;

    double Xi_r = Xi[i];

    double f_r_r1 = I1_r_r1[i];
    double f_r_r2 = I1_r_r2[i];
    double f_r_r1p = I1_r_r1p[i];
    double f_r_r2p = I1_r_r2p[i];

    for (size_t ll=0; ll<n_l2; ll++) {
      int ell2 = l2[ll];
      double wig = gsl_sf_coupling_3j(2*l, 2*l_prime, 2*ell2, 0, 0, 0);
      if (wig!=0) {
  double t1 = (2*ell2+1)*pow(wig,2);
  double f_r_r1_r1p = I2_r1_r1p[ll][i];
  double f_r_r2_r2p = I2_r2_r2p[ll][i];
  double f_r_r2_r1p = I2_r2_r1p[ll][i];
  double f_r_r1_r2p = I2_r1_r2p[ll][i];

  double t2 = pow(-1, l2[ll])*Xi_r*(f_r_r1_r1p*f_r_r2_r2p+f_r_r2_r1p*f_r_r1_r2p);
  double t3  = pow(-1, 0.5*(l+l_prime+ell2))*(f_r_r1*f_r_r1p*f_r_r2_r2p+f_r_r1*f_r_r2p*f_r_r2_r1p+f_r_r2*f_r_r1p*f_r_r1_r2p+f_r_r2*f_r_r2p*f_r_r1_r1p);
  sum += t1*(t2+t3);
      }
    }

    Int += rr[i]*rr[i]*sum*rr[i];
  }
  
  double fact = 4.*par::pi/Volume*(2*l+1)*(2*l_prime+1)*pow(-1, l+l_prime);
  
  return fact*Int*dlogr;
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::ThreePointCorrelation::zeta_covariance (const double Volume, const double nObjects, const std::vector<double> theta, const double r1, const double r2, const double deltaR, const std::vector<double> kk, const std::vector<double> Pk, const int norders, const double prec, const bool method, const int nExtractions, std::vector<double> mean, const int seed)
{
  (void)method;
  (void)nExtractions;
  (void)mean;
  (void)seed;

  vector<double> rr, Xi;
  wrapper::fftlog::transform_FFTlog(rr, Xi, 1, kk, Pk, 0);

  vector<vector<double>> zeta_l1l2_covariance(norders, vector<double>(norders, 0.));

  for (int i=0; i<norders; i++)
    for (int j=i; j<norders; j++) {
      zeta_l1l2_covariance[i][j] = zeta_multipoles_covariance(Volume, nObjects, i, j, r1, r2, r1, r2, deltaR, kk, Pk, rr, Xi, prec);
      zeta_l1l2_covariance[j][i] = zeta_l1l2_covariance[i][j];
    }

  const int ntheta = int(theta.size());

  vector<vector<double>> Pl_theta(ntheta, vector<double>(norders, 0));
  for (int i=0; i<ntheta; i++) {
    for (int j=0; j<norders; j++)
      Pl_theta[i][j] = legendre_polynomial (cos(theta[i]), j);
  }

  vector<vector<double>> zeta_covariance(ntheta, vector<double>(ntheta, 0));
  for (int i=0; i<ntheta; i++)
    for (int j=0; j<ntheta; j++)
      for (int l1=0; l1<norders; l1++)
  for (int l2=0; l2<norders; l2++)
    zeta_covariance[i][j] += zeta_l1l2_covariance[l1][l2]*Pl_theta[i][l1]*Pl_theta[j][l2];

  return zeta_covariance;
}


// =====================================================================================a


void cbl::cosmology::ThreePointCorrelation::xi_r_n (std::vector<double> &xi_n, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk)
{
  xi_n = wrapper::fftlog::transform_FFTlog (rr, 1, kk, Pk, nn, 0, par::pi, 1);
}


// =====================================================================================


void cbl::cosmology::ThreePointCorrelation::xi_r_n_pm (std::vector<double> &xi_n_p, std::vector<double> &xi_n_m, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk)
{
  vector<double> pk_p(Pk.size(), 0), pk_m(Pk.size(), 0);

  for (size_t i=0; i<Pk.size(); i++) {
    pk_p[i] = kk[i]*Pk[i];
    pk_m[i] = Pk[i]/kk[i];
  }

  xi_n_p = wrapper::fftlog::transform_FFTlog (rr, 1, kk, pk_p, nn, 0, par::pi, 1);
  xi_n_m = wrapper::fftlog::transform_FFTlog (rr, 1, kk, pk_m, nn, 0, par::pi, 1);
}


// =====================================================================================


void cbl::cosmology::ThreePointCorrelation::eff_l_l1 (std::vector<std::vector<double>> &eff, const std::vector<double> rr, const int l, const int l1, const std::vector<double> kk, const std::vector<double> Pk)
{
  double min_rr = Min(rr);
  double max_rr = Max(rr);
  vector<double> new_r = linear_bin_vector(rr.size(), min_rr, max_rr);
  eff.resize(rr.size());
  
  for (size_t i=0; i<rr.size(); i++)
  {
    vector<double> _pk(Pk.size(), 0);

    for (size_t j=0; j<Pk.size(); j++)
      _pk[j] = kk[j]*Pk[j]*jl(kk[j]*rr[i], l);

    eff[i] = wrapper::fftlog::transform_FFTlog (new_r, 1, kk, _pk, l1, 0, par::pi, 1);
  }
  
}


// =====================================================================================


void cbl::cosmology::ThreePointCorrelation::I_ELL_ell (std::vector<std::vector<double>> &II, const std::vector<double> rr, const int ll, const int LL, const std::vector<double> kk, const std::vector<double> Pk)
{
  II.resize(rr.size(), vector<double>(rr.size(), 0));
  double min_rr = Min(rr);
  double max_rr = Max(rr);
  vector<double> new_r = linear_bin_vector(rr.size(), min_rr, max_rr);

  for (int l1 = 0; l1<=LL+ll; l1++)
  {
    if ( (LL>= fabs(l1-ll)) && (LL <= l1+ll)) {
      double fact = pow(-1., l1+ll)*(2.*l1+1)*(2.*ll+1)*pow(gsl_sf_coupling_3j(2*l1, 2*ll, 2*LL, 0, 0, 0),2);

      if(fact!=0) {
  vector<vector<double>> eff;
  eff_l_l1 (eff, rr, ll, l1, kk, Pk);
  for (size_t r1=0; r1<rr.size(); r1++) {
    for (size_t r2=r1; r2<rr.size(); r2++) {

      glob::FuncGrid interp_r1_eff(new_r, eff[r1], "Spline");
      glob::FuncGrid interp_r2_eff(new_r, eff[r2], "Spline");

      auto integrand = [&] ( const double _r) {
        return interp_r1_eff(_r)*interp_r2_eff(_r)*_r;
      };
      II[r1][r2] += fact*wrapper::gsl::GSL_integrate_qag(integrand, min_rr, max_rr, 1.e-3); //Check the integral limits
      if(r1!=r2)
        II[r2][r1] += II[r1][r2];
      
    }
  }
      }
    }
  }
}


// =====================================================================================


void cbl::cosmology::ThreePointCorrelation::k_ell (std::vector<std::vector<double>> &KK, const std::vector<double> rr, const int ll, const std::vector<double> kk, const std::vector<double> Pk)
{

  vector<vector<double>> I1l, I3l, I5l;

  I_ELL_ell (I1l, rr, ll, 1, kk, Pk);
  I_ELL_ell (I3l, rr, ll, 3, kk, Pk);
  I_ELL_ell (I5l, rr, ll, 5, kk, Pk);

  KK.resize(rr.size(), vector<double>(rr.size(), 0));

  const double fact = 64./77175;
  for (size_t r1=0; r1<rr.size(); r1++)
    for (size_t r2=r1; r2<rr.size(); r2++) {
      KK[r1][r2] = fact*(9.*I1l[r1][r2]-14.*I3l[r1][r2]+5.*I5l[r1][r2]);

      if(r1!=r2)
  KK[r2][r1] = KK[r1][r2];
    }

}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_0_factor (const double b1, const double gamma, const double beta) const
{
  return pow(b1, 4)*beta*(2./3+38./45*beta+2./5*beta*beta+2./25*pow(beta, 3))+pow(b1, 3)*gamma*(1.+2./3*beta+1./9*beta*beta)+pow(b1, 3)*(34./21+94./63*beta+326./525*beta*beta+134./1225*pow(beta, 3));
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_0_factor_tidal (const double b1, const double gamma_t, const double beta) const
{
  return pow(b1, 3)*4.*gamma_t/225*(-75-50.*beta-7.*beta*beta);
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_1_factor (const double b1, const double beta) const
{
  return -pow(b1, 4)*beta*(1./3+3./5*beta+67./175*beta*beta+3./35*pow(beta, 3))-pow(b1, 3)*(1.+beta+37./75*beta*beta+17./175*pow(beta, 3));
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_2_factor (const double b1, const double gamma, const double beta) const
{
  return pow(b1, 4)*beta*beta*(16./45+16./35*beta+32./245*beta*beta)+pow(b1, 3)*(8./21+32./63*beta+144./245*beta*beta+4./45*gamma*beta*beta);  
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_2_factor_tidal (const double b1, const double gamma_t, const double beta) const
{
  return pow(b1, 3)*4.*gamma_t/63*(21.+14.*beta+beta*beta);
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_3_factor (const double b1, const double beta) const
{
  return -pow(b1, 4)*pow(beta, 3)*(8./175+8./315*beta)-pow(b1, 3)*beta*beta*(8./75+8./175*beta);
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_4_factor (const double b1, const double beta) const
{
  return pow(b1, 4)*128./11025*pow(beta, 4)+pow(b1, 3)*beta*beta*(-32./3675+32./8575*beta);
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_4_factor_tidal (const double b1, const double gamma_t, const double beta) const
{
  return pow(b1, 3)*32.*beta*beta*gamma_t/525;
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_ell_k_factor (const double b1, const double beta) const
{
  return pow(b1, 3)*(7.*beta*beta+3*pow(beta, 3));
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_precyclic_RSD (const double r1, const double r2, const double mu, const double b1, const double b2, const double bt, const double beta, std::vector<std::shared_ptr<cbl::glob::FuncGrid>> interp_xil, std::vector<std::shared_ptr<cbl::glob::FuncGrid>> interp_xil_m, std::vector<std::shared_ptr<cbl::glob::FuncGrid>> interp_xil_p, std::vector<std::shared_ptr<cbl::glob::FuncGrid2D>> interp_Kl, const bool use_k) const
{
  const double mu12 = mu;

  if (abs(mu12)>1)
    return 0;
    
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*mu12);
  const double a12 = acos(mu12);

  double a23 = 0.;
  double a13 = 0.;
  if (r1<=r2) {
    a23 = asin(sin(a12)*r1/r3);
    a13 = par::pi-a12-a23;
  }
  else {
    a13 = asin(sin(a12)*r2/r3);
    a23 = par::pi-a12-a13;
  }
  
  const double mu13 = cos(a13);
  const double mu23 = cos(a23);
  
  const double gamma = b2/b1;
  const double gamma_t = bt/b1;
  
  double fact = 0;
  
  // l = 0
  fact += (zeta_ell_0_factor(b1, gamma, beta)+zeta_ell_0_factor_tidal(b1, gamma_t, beta))*(interp_xil[0]->operator()(r1)*interp_xil[0]->operator()(r2));
  fact += (zeta_ell_0_factor(b1, gamma, beta)+zeta_ell_0_factor_tidal(b1, gamma_t, beta))*(interp_xil[0]->operator()(r1)*interp_xil[0]->operator()(r3));
  fact += (zeta_ell_0_factor(b1, gamma, beta)+zeta_ell_0_factor_tidal(b1, gamma_t, beta))*(interp_xil[0]->operator()(r2)*interp_xil[0]->operator()(r3));
  
  // l = 1
  fact += zeta_ell_1_factor(b1, beta)*(interp_xil_m[0]->operator()(r1)*interp_xil_p[0]->operator()(r2)+interp_xil_m[0]->operator()(r2)*interp_xil_p[0]->operator()(r1))*legendre_polynomial(mu12, 1);
  fact += zeta_ell_1_factor(b1, beta)*(interp_xil_m[0]->operator()(r1)*interp_xil_p[0]->operator()(r3)+interp_xil_m[0]->operator()(r3)*interp_xil_p[0]->operator()(r1))*legendre_polynomial(mu13, 1);
  fact += zeta_ell_1_factor(b1, beta)*(interp_xil_m[0]->operator()(r2)*interp_xil_p[0]->operator()(r3)+interp_xil_m[0]->operator()(r3)*interp_xil_p[0]->operator()(r2))*legendre_polynomial(mu23, 1);
  
  // l = 2
  fact += (zeta_ell_2_factor(b1, gamma, beta)+zeta_ell_2_factor_tidal(b1, gamma_t, beta))*(interp_xil[1]->operator()(r1)*interp_xil[1]->operator()(r2))*legendre_polynomial(mu12, 2);
  fact += (zeta_ell_2_factor(b1, gamma, beta)+zeta_ell_2_factor_tidal(b1, gamma_t, beta))*(interp_xil[1]->operator()(r1)*interp_xil[1]->operator()(r3))*legendre_polynomial(mu13, 2);
  fact += (zeta_ell_2_factor(b1, gamma, beta)+zeta_ell_2_factor_tidal(b1, gamma_t, beta))*(interp_xil[1]->operator()(r2)*interp_xil[1]->operator()(r3))*legendre_polynomial(mu23, 2);
  
  // l = 3
  fact += zeta_ell_3_factor(b1, beta)*(interp_xil_m[1]->operator()(r1)*interp_xil_p[1]->operator()(r2)+interp_xil_m[1]->operator()(r2)*interp_xil_p[1]->operator()(r1))*legendre_polynomial(mu12, 3);
  fact += zeta_ell_3_factor(b1, beta)*(interp_xil_m[1]->operator()(r1)*interp_xil_p[1]->operator()(r3)+interp_xil_m[1]->operator()(r3)*interp_xil_p[1]->operator()(r1))*legendre_polynomial(mu13, 3);
  fact += zeta_ell_3_factor(b1, beta)*(interp_xil_m[1]->operator()(r2)*interp_xil_p[1]->operator()(r3)+interp_xil_m[1]->operator()(r3)*interp_xil_p[1]->operator()(r2))*legendre_polynomial(mu23, 3);
  
  // l = 4
  fact += (zeta_ell_4_factor(b1, beta) + zeta_ell_4_factor_tidal(b1, gamma_t, beta))*(interp_xil[2]->operator()(r1)*interp_xil[2]->operator()(r2))*legendre_polynomial(mu12, 4);
  fact += (zeta_ell_4_factor(b1, beta) + zeta_ell_4_factor_tidal(b1, gamma_t, beta))*(interp_xil[2]->operator()(r1)*interp_xil[2]->operator()(r3))*legendre_polynomial(mu13, 4);
  fact += (zeta_ell_4_factor(b1, beta) + zeta_ell_4_factor_tidal(b1, gamma_t, beta))*(interp_xil[2]->operator()(r2)*interp_xil[2]->operator()(r3))*legendre_polynomial(mu23, 4);
  

  if (use_k) {
    for (int i=0; i<5; i++) {
      fact += zeta_ell_k_factor(b1, beta)*interp_Kl[i]->operator()(r1, r2)*legendre_polynomial(mu12, i);
      fact += zeta_ell_k_factor(b1, beta)*interp_Kl[i]->operator()(r1, r3)*legendre_polynomial(mu13, i);
      fact += zeta_ell_k_factor(b1, beta)*interp_Kl[i]->operator()(r2, r3)*legendre_polynomial(mu23, i);
      }
    }

  return fact;
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::zeta_expansion_RSD (const double r1, const double r2, const double b1, const double b2, const double bt, const double beta, std::vector<double> &rr, std::vector<std::vector<double>> &xil, std::vector<std::vector<double>> &xil_m, std::vector<std::vector<double>> &xil_p, std::vector<std::vector<std::vector<double>>> &Kl, const int norders, const double prec, const bool use_k) const
{   
  std::vector<std::shared_ptr<glob::FuncGrid>> interp_xil;
  for (const auto& x : xil) 
    interp_xil.push_back(std::make_shared<glob::FuncGrid>(rr, x, "Spline"));
  
  std::vector<std::shared_ptr<glob::FuncGrid>> interp_xil_m;
  for (const auto& x : xil_m) 
    interp_xil_m.push_back(std::make_shared<glob::FuncGrid>(rr, x, "Spline"));
  
  std::vector<std::shared_ptr<glob::FuncGrid>> interp_xil_p;
  for (const auto& x : xil_p) 
    interp_xil_p.push_back(std::make_shared<glob::FuncGrid>(rr, x, "Spline"));
  
  std::vector<std::shared_ptr<glob::FuncGrid2D>> interp_Kl;
  for (const auto& x : Kl) 
    interp_Kl.push_back(std::make_shared<glob::FuncGrid2D>(rr, rr, x, "Linear"));

  vector<double> zeta_r1_r2(norders, 0);

  for (int i=0; i<norders; i++) {
    auto integrand = [&] (const double mu12) {
      return ThreePointCorrelation::zeta_precyclic_RSD(r1, r2, mu12, b1, b2, bt, beta, interp_xil, interp_xil_m, interp_xil_p, interp_Kl, use_k)*legendre_polynomial (mu12, i);
    };
    zeta_r1_r2[i] = 0.5*(2*i+1)*wrapper::gsl::GSL_integrate_qag(integrand, -1, 1, prec);
  }

  return zeta_r1_r2;
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &xi_matter_m1, std::vector<double> &xi_matter_p1, std::vector<double> &xi_matter_2, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders, const double prec, const bool use_k) const
{
  if (rr.size()==0) {
    rr = linear_bin_vector(200, 1., 300.);
    integrals_zeta_Slepian(xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, rr, kk, Pk_matter);
  }

  const double mu = cos(theta);
  vector<double> z_r1_r2 = ThreePointCorrelation::zeta_expansion_RSD(r1, r2, 1., 0., 0., 0., rr, xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, norders, prec, use_k);

  double zeta_r1_r2_theta = 0.;
  for (size_t i=0; i<z_r1_r2.size(); i++)
    zeta_r1_r2_theta += z_r1_r2[i]*legendre_polynomial(mu, i);

  return zeta_r1_r2_theta;
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::zeta_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_matter;
  vector<double> zH(ntheta, 0);

  if (model=="Slepian") {
    vector<double> xi_matter_m1, xi_matter_p1, xi_matter_2;
    for (int i=0; i<ntheta; i++)
      zH[i] = b1*b1*b1*zeta_DM_Slepian(r1, r2, theta[i], rr, xi_matter, xi_matter_m1, xi_matter_p1, xi_matter_2, kk, Pk_matter, 9, 1.e-3)+b1*b1*b2*denominator_Q(r1, r2, theta[i], rr, xi_matter);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<ntheta; i++)
      zH[i] = b1*b1*b1*zeta_DM_BarrigaGatzanaga(r1, r2, theta[i], rr, xi_matter, Phi, kk, Pk_matter)+b1*b1*b2*denominator_Q(r1, r2, theta[i], rr, xi_matter);
  }
  else
    ErrorCBL("the chosen model is not implemented!", "z_halo", "3PCF.cpp");

  return zH;
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_DM_Slepian_RSD (const double r1, const double r2, const double theta, const double beta, std::vector<double> &rr, std::vector<std::vector<double>> &xil, std::vector<std::vector<double>> &xil_m, std::vector<std::vector<double>> &xil_p, std::vector<std::vector<std::vector<double>>> &Kl, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders, const double prec, const bool use_k)
{
  if (rr.size()==0) {
    rr = linear_bin_vector(200, 1., 300.);
    for (int i=0; i<3; i++) {
      double ll = 2*i;
      xi_r_n(xil[i], rr, ll, kk, Pk_matter);
      }
    for (int i=0; i<2; i++) {
      double ll = 2*i+1;
      xi_r_n_pm(xil_p[i], xil_m[i], rr, ll, kk, Pk_matter);
      }
    for (int i=0; i<5; i++)
      //k_ell(Kl[i], rr, i, kk, Pk_matter);
      Kl[i].resize(rr.size(), std::vector<double>(rr.size(), 0.));
  }

  const double mu = cos(theta);
  vector<double> z_r1_r2 = ThreePointCorrelation::zeta_expansion_RSD(r1, r2, 1, 0, 0, beta, rr, xil, xil_m, xil_p, Kl, norders, prec, use_k);

  double zeta_r1_r2_theta = 0.;
  for (size_t i=0; i<z_r1_r2.size(); i++)
    zeta_r1_r2_theta += z_r1_r2[i]*legendre_polynomial(mu, i);

  return zeta_r1_r2_theta;
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::zeta_halo_RSD (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double bt, const double beta, const string model, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders, const bool use_k)
{
  const int ntheta = theta.size();
  vector<double> rr;
  vector<double> zH(ntheta, 0);

  if (model=="Slepian") {
    std::vector<std::vector<double>> xil(3);
    std::vector<std::vector<double>> xil_m(2);
    std::vector<std::vector<double>> xil_p(2);
    std::vector<std::vector<std::vector<double>>> Kl(5);
    const double prec = 1.e-3;

    for (int i=0; i<ntheta; i++)
      zH[i] = b1*b1*b1*zeta_DM_Slepian_RSD(r1, r2, theta[i], beta, rr, xil, xil_m, xil_p, Kl, kk, Pk_matter, norders, prec, use_k)+b1*b1*b2*denominator_Q(r1, r2, theta[i], rr, xil[0])+bt*0;
  }
  else
    ErrorCBL("the chosen model is not implemented!", "z_halo_RSD", "3PCF.cpp");

  return zH;
}


// =====================================================================================


double cbl::cosmology::ThreePointCorrelation::zeta_RSD (const double r1, const double r2, const double theta, const double b1, const double b2, const double bt, const double beta, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders, const double prec, const bool use_k)
{
  double zeta_g = 0;
  std::vector<double> rr = linear_bin_vector(200, 1., 300.);
  std::vector<std::vector<double>> xil(3);
  std::vector<std::vector<double>> xil_m(2);
  std::vector<std::vector<double>> xil_p(2);
  std::vector<std::vector<std::vector<double>>> Kl(5);

   for (int i=0; i<3; i++) {
    double ll = 2*i;
    xi_r_n(xil[i], rr, ll, kk, Pk_matter);
    }
  for (int i=0; i<2; i++) {
    double ll = 2*i+1;
    xi_r_n_pm(xil_p[i], xil_m[i], rr, ll, kk, Pk_matter);
    }
  for (int i=0; i<5; i++)
    //k_ell(Kl[i], rr, i, kk, Pk_matter);                                                                                                                                                                   
    Kl[i].resize(rr.size(), std::vector<double>(rr.size(), 0.));

  std::vector<double> z_ell = ThreePointCorrelation::zeta_expansion_RSD(r1, r2, b1, b2, bt, beta, rr, xil, xil_m, xil_p, Kl, norders, prec, use_k);

  double mu = cos(theta);
  for (size_t j=0; j<z_ell.size(); j++)
    zeta_g += z_ell[j]*legendre_polynomial(mu, j);

  return zeta_g;
}


// =====================================================================================


std::vector<double> cbl::cosmology::ThreePointCorrelation::zeta_RSD (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double bt, const double beta, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders, const double prec, const bool use_k)
{
  const int ntheta = theta.size();
  std::vector<double> zeta_g(ntheta, 0.);

  for (int i=0; i<ntheta; i++)
    zeta_g[i] = zeta_RSD(r1, r2, theta[i], b1, b2, bt, beta, kk, Pk_matter, norders, prec, use_k);

  return zeta_g;
}

