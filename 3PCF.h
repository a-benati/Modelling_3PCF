/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *******************************************************************/

/**
 *  @file Headers/3PCF.h
 *
 *  @brief The class ThreePointCorrelation
 *
 *  This file defines the interface of the class
 *  ThreePointCorrelation, used to model three-point statistics
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __THREEPOINTCORRELATION__
#define __THREEPOINTCORRELATION__

#include "PkXi.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class ThreePointCorrelation 3PCF.h "Headers/3PCF.h"
     *
     *  @brief The class ThreePointCorrelation
     *
     *  This class is used to handle objects of type <EM>
     *  ThreePointCorrelation </EM>. It is used to model three-point
     *  statistics
     */
    class ThreePointCorrelation {

    private:
      
      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      ThreePointCorrelation () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      ThreePointCorrelation (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~ThreePointCorrelation () = default;

      ///@}

      
      /**
       *  @name Functions to get the private members of the class
       */
      ///@{
      
      /**
       *  @brief Get the private member m_cosmology
       *
       *  @return the cosmological model
       */
      std::shared_ptr<cosmology::Cosmology> cosmology ()
      { return move(m_cosmology); }

      ///@}
      

      /**
       *  @name Functions to estimate the three-point correlation
       *  model
       */
      ///@{

      /**
       *  @brief the normalization factor for reduced three-point
       *  correlation function
       *
       *  this function computes the normalization factor for reduced
       *  three-point correlation function:
       *
       *  \f[ \xi(r_1)\cdot\xi(r_2) + \xi(r_2)\cdot\xi(r_3) +
       *  \xi(r_3)\cdot\xi(r_1) \f]
       *
       *  with \f$ r_3 = \sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos(\theta)}
       *  \f$
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle between r1 and r2
       *
       *  @param rr vector containing the scales at which the
       *  two-point correlation function is computed
       *  
       *  @param xi_matter vector containing the dark matter two-point
       *  correlation function values, estimated at the scales given
       *  in rr
       *
       *  @return the normalization factor for reduced three-point
       *  correlation function
       */
      double denominator_Q (const double r1, const double r2, const double theta, const std::vector<double> rr, const std::vector<double> xi_matter) const;

      /**
       *  @brief integral functions for the three-point correlation
       *  model
       *
       *  this function computes and store functons used to model the
       *  three-point correlation model; specifically, it implements
       *  Eq. 21, in polar coordinates, of Bel et al. 2015, MNRAS,
       *  453, 259):
       *
       *  \f[ \xi_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty \mathrm{d}
       *  k\, k^2 P_{DM}(k) j_0(k r), \\\
       *
       *  \Phi(r) = \frac{1}{2\pi^2}\int_0^\infty \mathrm{d} k\,
       *  P_{DM}(k) W^2(kr) j_0(k r), \\ \f]
       *
       *  where \f$j_0(k r)=\sin(k r)/(kr)\f$ is the l=0 spherical
       *  Bessel function, and \f$W(kr)\f$ is the top-hat window
       *  function computed by cbl::TopHat_WF
       *
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point correlation function values
       *
       *  @param [out] Phi vector containing the \f$ \Phi(r)\f$
       *  values, estimated at the scales given in rr
       *
       *  @param [in] rr vector or scales at which the dark matter
       *  two-point correlation function (xi_matter) will be computed
       *
       *  @param [in] kk vector of the wave vector modules at which the
       *  power spectrum is computed
       *
       *  @param [in] Pk_matter vector of containing the dark matter power
       *  spectrum values, estimated at the wave vector modules given
       *  in kk
       *
       *  @param prec the integral precision
       */
      void integrals_Q_nonLocal (std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_matter, const double prec) const;

      /**
       *  @brief function to compute non-local contribution to
       *  three-point correlation function; specifically, it
       *  implements Eq. 20 of Bel el et al. 2015, MNRAS, 453, 259:
       *
       *  \f[ \Gamma_{123} = \left[
       *  \xi(r_1)+3\frac{\Phi^\prime(r_1)}{r1}\right] \left[
       *  \xi(r_2)+3\frac{\Phi^\prime(r_2)}{r_2}\right]P_2(\cos\theta)
       *  \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, \f$P_2\f$ is the second Legandre polynomial
       *  computed by cbl::legendre_polynomial, and \f$\xi(r),
       *  \Phi(r)\f$ are the integrals of the power spectrum computed
       *  by cbl::cosmology::Cosmology::integrals_Q_nonLocal
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle betwee r1 and r2
       *
       *  @param xi vector containing the value of xi at r1, r2
       *
       *  @param dPhi vector containing the value of the derivative
       *  of Phi at r1, r2
       *
       *  @return the value of the \f$\Gamma_{123}\f$
       */
      double Gamma_3PCF (const double r1, const double r2, const double theta, const std::vector<double> xi, const std::vector<double> dPhi) const;

      /**
       *  @brief the non-local contribution to the reduced dark
       *  matter three-point correlation function
       *
       *  this function computes the non-local contribution to
       *  three-point correlation function; specifically, it
       *  implements Eq. 22 of Bel el et al. 2015, MNRAS, 453, 259:
       *
       *  \f[ Q_{non-local}(r_1, r_2, \theta) = \frac{2}{3} \left(
       *  \frac{\Gamma_{123} + \Gamma_{312} + \Gamma_{231}}
       *  {\xi(r_1)\cdot\xi(r_2) + \xi(r_2)\cdot\xi(r_3) +
       *  \xi(r_3)\cdot\xi(r_1)}-1 \right) \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, and \f$\xi(r), \Phi(r)\f$ are the integrals of the
       *  power spectrum computed by
       *  cbl::cosmology::Cosmology::integrals_Q_nonLocal
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle betwee r1 and r2
       *
       *  @param [out] rr vector or scales at which the dark matter
       *  two-point correlation function is computed
       *
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point correlation function values
       *
       *  @param [out] Phi vector containing the \f$ \Phi(r)\f$
       *  values, estimated at the scales given in rr
       *
       *  @param [in] kk vector of the wave vector modules at which
       *  the power spectrum is computed
       *
       *  @param [in] Pk_matter vector of containing the dark matter
       *  power spectrum values, estimated at the wave vector modules
       *  given in kk
       *
       *  @return the value of non-local Q
       */
      double Q_nonLocal (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief all the non-local contribution terms of the reduced
       *  dark matter three-point correlation function
       *  
       *  this function computes all the the non-local contribution
       *  terms of the reduced three-point correlation function,
       *  computed by cbl::cosmology::Cosmology::Q_nonLocal
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle betwee r1 and r2
       *
       *  @param kk vector of the wave vector modules at which the
       *  power spectrum is computed
       *
       *  @param Pk_matter vector of containing the dark matter power
       *  spectrum values, estimated at the wave vector modules given
       *  in kk
       *
       *  @return vector containing the DM reduced three-point
       *  correlation function
       */
      std::vector<double> Q_nonLocal (const double r1, const double r2, const std::vector<double> theta, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief integrals used to compute the Slepian et al. 2015
       *  three-point correlation function model
       *
       *  this function computes the integrals used to
       *  model the three-point correlation as described in Slepian
       *  et. al 2015:
       *
       *  \f[ \xi_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty
       *  \mathrm{d} k k^2 P_{DM}(k) j_0(k r), \\\
       *  \xi^{[1\pm]}_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty
       *  \mathrm{d} k k^2 P_{DM}(k) k^{\pm 1} j_1(k r), \\
       *  \xi^{[2]}_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty
       *  \mathrm{d} k k^2 P_{DM}(k) j_2(k r) . \f]
       *
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] xi_matter_m1 vector containing
       *  \f$\xi^{[1-]}_{DM}(r)\f$
       *
       *  @param [out] xi_matter_p1 vector containing
       *  \f$\xi^{[1+]}_{DM}(r)\f$
       *
       *  @param [out] xi_matter_2 vector containing
       *  \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @param [in] rr vector or scales
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_matter the dark matter power spectrum
       */
      void integrals_zeta_Slepian (std::vector<double> &xi_matter, std::vector<double> &xi_matter_m1, std::vector<double> &xi_matter_p1, std::vector<double> &xi_matter_2, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the dark matter three-point correlation function
       *  model by Slepian et al. 2015
       *
       *  this function computes \f$\zeta_{DM} (r_1, r_2, \hat{r_1}
       *  \cdot \hat{r_2})\f$, as described in Slepian et al. 2015:
       *
       *  \f[ \zeta_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = 
       *  \sum_l \zeta_l(r_1, r_2) P_l(\hat{r_1} \cdot \hat{r_2}) .\f]
       *
       *  The coefficients of the expansion are computed by
       *  cbl::cosmology::Cosmology::zeta_expansion_RSD (using \f$b_1 = 1\f$, \f$b_2 = 0\f$, \f$b_t = 0\f$ and \f$\beta = 0\f$)
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] xi_matter_m1 vector containing
       *  \f$\xi^{[1-]}_{DM}(r)\f$
       *
       *  @param [out] xi_matter_p1 vector containing
       *  \f$\xi^{[1+]}_{DM}(r)\f$
       *
       *  @param [out] xi_matter_2 vector containing
       *  \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_matter the dark matter power spectrum
       *
       *  @param [in] norders the maximum number of orders
       *
       *  @param [in] prec the integral precision
       *
       *  @return the connected dark matter three-point correlation
       *  function
       */
      double zeta_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &xi_matter_m1, std::vector<double> &xi_matter_p1, std::vector<double> &xi_matter_2, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders=9, const double prec=1.e-3, const bool use_k=false) const;

      /**
       *  @brief the dark matter reduced three-point correlation
       *  function model by Slepian et al. 2015
       *
       *  this function computes \f$Q_{DM} (r_1, r_2, \hat{r_1} \cdot
       *  \hat{r_2})\f$ as described in Slepian et al. 2015:
       *
       *  \f[ Q_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{\zeta_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2})}
       *  {\left(
       *  \xi(r_1)\xi(r_2)+\xi(r_2)\xi(r_3)+\xi(r_3)\xi(r_1)\right)}
       *  \f]
       *
       *  see cbl::cosmology::Cosmology::zeta_DM_Slepian for
       *  more details
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point orrelation function
       *
       *  @param [out] xi_matter_m1 vector containing
       *  \f$\xi^{[1-]}_{DM}(r)\f$
       *
       *  @param [out] xi_matter_p1 vector containing
       *  \f$\xi^{[1+]}_{DM}(r)\f$
       *
       *  @param [out] xi_matter_2 vector containing
       *  \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @param [out] kk vector of the wave vector modules
       *
       *  @param [out] Pk_matter the dark matter power spectrum
       *
       *  @param [in] norders the maximum numbers of orders
       *
       *  @param [in] prec the integral precision                                                                                                                                                       
       *
       *  @param [in] use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return the dark matter reduced three-point correlation
       *  function
       */
      double Q_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &xi_matter_m1, std::vector<double> &xi_matter_p1, std::vector<double> &xi_matter_2, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders=9, const double prec=1.e-3) const;

      /**
       *  @brief integrals used to compute the Barriga & Gatzanaga
       *  al. 2002 three-point correlation function model
       * 
       *  this function computes the integrals used to model the
       *  three-point correlation as described in Barriga & Gatzanaga
       *  2002:
       *
       *  \f[ \xi_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty \mathrm{d}
       *  k\, k^2 P_{DM}(k) j_0(k r), \f] \f[ \Phi(r) =
       *  \frac{1}{2\pi^2}\int_0^\infty \mathrm{d} k\, P_{DM}(k)
       *  j_0(k r). \f]
       *
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] Phi vector containing \f$ \Phi(r)\f$
       *
       *  @param [in] rr vector or scales
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_matter the dark matter power spectrum
       */
      void integrals_zeta_BarrigaGatzanaga (std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the single term of the dark matter three-point
       *  correlation function model by Barriga & Gatzanaga et
       *  al. 2002
       *
       *  this function computes the single term of the dark matter
       *  three-point correlation function, following Barriga &
       *  Gatzanaga et al.  2002:
       *
       *  \f[ f(r_1, r_2) = \frac{10}{7}\xi(r_1) \xi(r_2)+\frac{4}{7}
       *  \left\{ -3 \frac{\Phi^\prime(r_1) \Phi^\prime(r_2)}{r_1
       *  r_2} -\frac{\xi(r_1) \Phi^\prime(r_2)}{r_2}-\frac{\xi(r_2)
       *  \Phi^\prime(r_1)}{r_1} +\mu^2\left[
       *  \xi(r_1)+3\frac{\Phi^\prime(r_1)}{r1}\right]\left[
       *  \xi(r_2)+3\frac{\Phi^\prime(r_2)}{r_3}\right] \right\}
       *  -\mu\left[ \xi^\prime(r_1)\Phi^\prime(r_2) +
       *  \xi^\prime(r_2)\Phi^\prime(r_1)\right] \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, and \f$\xi(r), \Phi(r)\f$ are the integrals of the
       *  power spectrum computed by
       *  cbl::cosmology::Cosmology::integrals_zeta_BarrigaGatzanaga
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle betwee r1 and r2
       *
       *  @param xi vector containing the value of xi at r1, r2
       *
       *  @param dxi vector containing the value of the derivative of
       *  xi at r1, r2
       *
       *  @param dPhi vector containing the value of the derivative
       *  of Phi at r1, r2
       *
       *  @return the dark matter reduced three-point correlation
       *  function
       */
      double zeta_single_BarrigaGatzanaga (const double r1, const double r2, const double theta, const std::vector<double> xi, const std::vector<double> dxi, const std::vector<double> dPhi) const;

      /**
       *  @brief the dark matter three-point correlation function
       *  model by Barriga & Gatzanaga et al. 2002
       *
       *  this functions computes the dark matter three-point
       *  correlation function model by Barriga & Gatzanaga et al
       *  2002:
       *
       *  \f[ f(r_1, r_2) = \frac{10}{7}\xi(r_1) \xi(r_2)+\frac{4}{7}
       *  \left\{ -3 \frac{\Phi^\prime(r_1) \Phi^\prime(r_2)}{r_1
       *  r_2} -\frac{\xi(r1) \Phi^\prime(r_2)}{r_2}-\frac{\xi(r2)
       *  \Phi^\prime(r_1)}{r_1} +\mu^2\left[
       *  \xi(r_1)+3\frac{\Phi^\prime(r_1)}{r1}\right]\left[
       *  \xi(r_2)+3\frac{\Phi^\prime(r_2)}{r_3}\right] \right\}
       *  -\mu\left[ \xi^\prime(r_1)\Phi^\prime(r_2) +
       *  \xi^\prime(r_2)\Phi^\prime(r_1)\right] +
       *  \mathrm{permutations} \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, and \f$\xi(r), \Phi(r)\f$ are the integrals of the
       *  power spectrum computed by
       *  cbl::cosmology::Cosmology::integrals_zeta_BarrigaGatzanaga.
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] Phi vector containing \f$ \Phi(r)\f$
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_matter the dark matter power spectrum
       *
       *  @return the dark matter three-point correlation function
       */
      double zeta_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the dark matter reduced three-point correlation
       *  function model by Barriga & Gatzanaga et al. 2002
       *
       *  this functions computes \f$Q_{DM} (r_1, r_2, \hat{r_1}
       *  \cdot \hat{r_2})\f$, as described in Barriga & Gatzanaga et
       *  al. 2002:
       *
       *  \f[ Q_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{\zeta_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2})}
       *  {\left(
       *  \xi(r_1)\xi(r_2)+\xi(r_2)\xi(r_3)+\xi(r_3)\xi(r_1)\right)}
       *  \f]
       *
       *  see cbl::cosmology::Cosmology::zeta_DM_BarrigaGatzanaga
       *  for more details
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       * 
       *  @param [out] rr vector or scales
       * 
       *  @param [out] xi_matter vector containing the dark matter
       *  two-point correlation function
       * 
       *  @param [out] Phi vector containing \f$ \Phi(r)\f$
       * 
       *  @param [in] kk vector of the wave vector modules
       * 
       *  @param [in] Pk_matter the dark matter power spectrum
       *
       *  @return the dark matter reduced three-point correlation
       *  function
       */
      double Q_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_matter, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the dark matter three-point correlation function
       *
       *  this function computes the dark matter three-point
       *  correlation function with either the Slepian et al 2015 or
       *  the Barriga & Gatzagnaga 2002 model
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @return vector containing the dark matter three-point
       *  correlation function
       */
      std::vector<double> zeta_DM (const double r1, const double r2, const std::vector<double> theta, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the dark matter reduced three-point correlation
       *  function
       *
       *  this function computes the dark matter reduced three-point
       *  reduced correlation function with either the Slepian et al.
       *  2015 or the Barriga & Gatzagnaga 2002 model
       *
       *  @param r1 the first side of the triangle 
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @return vector containing the dark matter reduced
       *  three-point correlation function
       */
      std::vector<double> Q_DM (const double r1, const double r2, const std::vector<double> theta, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the local-bias model of the three-point correlation
       *  function of dark matter haloes in real space
       *
       *  this function computes the three-point correlation function
       *  of dark matter haloes  in real space with either the Slepian et al.  2015
       *  or the Barriga & Gatzagnaga 2002 model, as follows:
       *
       *  \f[ \zeta_h (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = b_1^3
       *  \zeta_{DM}(r_1, r_2, \hat{r_1} \cdot \hat{r_2}) + b_1^2 b_2
       *  \left[ \xi(r_1)\cdot\xi(r_2) + \xi(r_2)\cdot\xi(r_3) +
       *  \xi(r_3)\cdot\xi(r_1) \right] \f]
       *
       *  with \f$r_3 = \sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos(\theta)}\f$
       *  and \f$b_1, b_2\f$ the linear and non-linear halo bias,
       *  respectively; \f$\zeta_{DM}\f$ is compute by
       *  cbl::cosmology::Cosmology::zeta_DM
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *  
       *  @param theta vector containing angles between r1 and r2, in
       *  radians 
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @return vector containing the three-point correlation
       *  function of dark matter haloes in real space
       */
      std::vector<double> zeta_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the local-bias model of the reduced three-point
       *  correlation function of dark matter haloes
       *
       *  this function computes the reduced three-point correlation
       *  function of dark matter haloes with either the Slepian et
       *  al.  2015 or the Barriga & Gatzagnaga 2002 model, as
       *  follows:
       *
       *  \f[ Q_h (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{Q_{DM}(r_1, r_2, \hat{r_1} \cdot
       *  \hat{r_2})}{b_1}+\frac{b_2}{b_1^2}\f]
       *
       *  \f$Q_{DM}\f$ is compute by
       *  cbl::cosmology::Cosmology::Q_DM
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @return vector containing the reduced three-point
       *  correlation function of dark matter haloes
       */
      std::vector<double> Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the non-local-bias model of the three-point
       *  correlation function of dark matter haloes
       *
       *  this function computes the reduced three-point correlation
       *  function of dark matter haloes, with non-local bias
       *  corrections, with either the Slepian et al. 2015 or the
       *  Barriga & Gatzagnaga 2002 model, as follows:
       *
       *  \f[ Q_h (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{Q_{DM}(r_1, r_2, \hat{r_1} \cdot \hat{r_2})}{b_1}+
       *  \frac{b_2}{b_1^2}+\frac{g_2}{b_1}Q_{non-local} \f]
       *
       *  \f$Q_{DM}\f$ is compute by
       *  cbl::cosmology::Cosmology::Q_DM and \f$Q_{non-local}\f$
       *  is the non-local contirbuion term, computed by
       *  cbl::cosmology::Cosmology::Q_nonLocal
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param g2 the non-local bias
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules 
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @return vector containing the reduced three-point
       *  correlation function of dark matter haloes
       */
      std::vector<double> Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double g2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the dark matter equilateral three-point correlation
       *  function
       *
       *  this function computes the dark matter equilateral
       *  three-point correlation function with either the Slepian et
       *  al 2015 or the Barriga & Gatzagnaga 2002 model
       *
       *  @param rr vector of sides
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @return vector containing the dark matter three-point
       *  correlation function
       */
      std::vector<double> zeta_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the dark matter equilateral reduced three-point
       *  correlation function
       *
       *  this function computes the dark matter equilateral reduced
       *  three-point correlation function with either the Slepian et
       *  al 2015 or the Barriga & Gatzagnaga 2002 model
       *
       *  @param rr vector of sides
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @return vector containing the dark matter three-point
       *  correlation function
       */
      std::vector<double> Q_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter) const;

      /**
       *  @brief the dark matter three-point correlation function 
       *  multipoles covariance model, by Slepian et al. 2015
       *
       *  \f[
       *  C_{{\rm GRF}, ll'}(r_{1},r_{2};r_{1}',r_{2}')=\frac{4\pi}{V}(2l+1)(2l'+1)(-1)^{l+l'} \times\int r^{2}dr\sum_{l_{2}}(2l_{2}+1)
       *  \left(\begin{array}{ccc} 
       *  l & l' & l_{2}\\
       *  0 & 0 & 0
       *   \end{array}\right)^2\nonumber\\
       *   \times\bigg\{(-1)^{l_2}\xi(r)\bigg[f_{l_{2}ll'}(r;r_{1},r_{1}')f_{l_{2}ll'}(r;r_{2},r_{2}')
       *   +f_{l_{2}ll'}(r;r_{2},r_{1}')f_{l_{2}ll'}(r;r_{1},r_{2}')\bigg]+(-1)^{(l+l'+l_{2})/2}
       *   \times\bigg[f_{ll}(r;r_{1})f_{l'l'}(r;r_{1}')f_{l_{2}ll'}(r;r_{2},r_{2}')
       *   +f_{ll}(r;r_{1})f_{l'l'}(r;r_{2}')f_{l_{2}ll'}(r;r_{2},r_{1}')
       *   +f_{ll}(r;r_{2})f_{l'l'}(r;r_{1}')f_{l_{2}ll'}(r;r_{1},r_{2}')
       *   +f_{ll}(r;r_{2})f_{l'l'}(r;r_{2}')f_{l_{2}ll'}(r;r_{1},r_{1}')\bigg]\bigg\}
       *  \f]
       *
       *  with:
       *  \f[
       *   f_{ll}(r;r_{1})=\int\frac{k^{2}dk}{2\pi^{2}}\left[P(k)+\frac{1}{n}\right]j_{l}(kr_{1})j_{l}(kr)
       *  \f]
       *   and:
       *   \f[
       *   f_{l_{2}ll'}(r;r_{1},r_{1}')=\int\frac{k^{2}dk}{2\pi^{2}} \left[P(k)+\frac{1}{n}\right] j_{l}(kr_{1})j_{l'}(kr_{1}')j_{l_{2}}(kr),
       *   \f]
       *   where \f$V\f$ is the effective survey volume, and  \f$n\f$ is the effective number density of the survey.
       *
       *  @param Volume the volume
       *
       *  @param nObjects the number of objects
       *
       *  @param l the order l of the multipoles expansion
       *
       *  @param l_prime the order \f$l'\f$ of the multipoles
       *  expansion
       *
       *  @param r1 the scale \f$r_1\f$
       *  
       *  @param r2 the scale \f$r_2\f$
       *  
       *  @param r1_prime the scale \f$r_1'\f$
       *
       *  @param r2_prime the scale \f$r_2'\f$
       *
       *  @param deltaR the bin size, if non-positive, no bin average 
       *  is computed 
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk the pdark matter ower spectrum
       *
       *  @param rr vector of scales
       *
       *  @param Xi vector containing the two-point correlation
       *  function
       *
       *  @param prec the integral precision
       *
       *  @return the covariance of the multipole expansion of 
       *  the dark matter three-point correlation function
       */
      double zeta_multipoles_covariance (const double Volume, const double nObjects, const int l, const int l_prime, const double r1, const double r2, const double r1_prime, const double r2_prime, const double deltaR, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> rr, const std::vector<double> Xi, const double prec=1.e-3);

      /**
       *  @brief the dark matter three-point correlation function 
       *  covariance model
       *
       *  this function computes the dark matter three-point
       *  correlation function covariance model, by Slepian et
       *  al. 2015, as a function of \f$\theta = r_1 \cdot r_2\f$:
       *
       *  \f[ C(r_1, r_2, \vec{r_1}\cdot\vec{r_2} \equiv \cos(\theta)
       *  = \sum_{l=0}^{l=max_l} \sum_{l^'=0}^{l^'=max_l} C_{l,
       *  l^'}(r_1, r_2) P_l(\cos(\theta) P_{l^'}(\cos(theta) \f]
       *
       *  where \f$C_{l, l^'}(r_1, r_2)\f$ is computed by
       *  cbl::cosmology::Cosmology::zeta_multipoles_covariance
       *
       *  @param Volume the volume
       *
       *  @param nObjects the number of objects
       *
       *  @param theta vector of angles at which the covariance is
       *  computed
       *
       *  @param r1 the scale \f$r_1\f$
       *  
       *  @param r2 the scale \f$r_2\f$
       *
       *  @param deltaR the bin size, if non-positive, no bin average 
       *  is computed 
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk the pdark matter ower spectrum
       *
       *  @param norders the maximum number of orders of multipoles
       *  of the three point correlation function expansion
       *
       *  @param prec the integral precision
       *
       *  @param method false \f$\rightarrow\f$ apply method 1; true
       *  \f$\rightarrow\f$ apply method 2
       *
       *  @param nExtractions the number of mock extraction
       *  from zeta multipoles coefficient covariance matrix
       *
       *  @param mean vector containing the mean values
       *
       *  @param seed random number generator seed
       *
       *  @return the covariance of the dark matter 
       *   three-point correlation function
       */
      std::vector<std::vector<double>> zeta_covariance (const double Volume, const double nObjects, const std::vector<double> theta, const double r1, const double r2, const double deltaR, const std::vector<double> kk, const std::vector<double> Pk, const int norders=10, const double prec=1.e-3, const bool method=false, const int nExtractions=10000, const std::vector<double> mean={}, const int seed=543);

      /**
       * @brief compute the  power spectrum integral transform
       *
       * this function computes the power spectrum integral transform:
       *
       * \f[ 
       *   \xi^{[n]} (r) = \int \frac{k^2\mathrm{d}k}{2\pi^2} P(k) j_n(kr).
       * \f]
       *
       * where n is the order of the transform.
       *
       * @param xi_n the power spectrum transform \f$\xi^{[n]} (r)\f$
       *
       * @param rr vector of scales
       *
       * @param nn the order of the transform
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       */
      void xi_r_n (std::vector<double> &xi_n, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the  power spectrum integral transform
       *
       * this function computes the power spectrum integral transform:
       *
       * \f[ 
       *   \xi^{[n\pm]} (r) = \int \frac{k^2\mathrm{d}k}{2\pi^2} k^{\pm1} P(k) j_n(kr)
       * \f]
       *
       * where n is the order of the transform.
       *
       * @param xi_n_p the power spectrum transform \f$\xi^{[n+]} (r)\f$
       * 
       * @param xi_n_m the power spectrum transform \f$\xi^{[n-]} (r)\f$
       *
       * @param rr vector of scales
       *
       * @param nn the order of the transform
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       */
      void xi_r_n_pm (std::vector<double> &xi_n_p, std::vector<double> &xi_n_m, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the  power spectrum integral transform
       *
       * this function computes the power spectrum integral transform:
       *
       * \f[ 
       *    f_{l, l_1} (r_i;r)= \int \frac{k^2\mathrm{d}k}{2\pi^2} j_l(kr_i) j_{l_1} (kr) k P(k)
       * \f]
       *
       * where \f$l, l_1\f$ are the orders of the transform.
       *
       * @param eff the power spectrum transform \f$ f_{l, l_1} (r_i;r) \f$
       * 
       * @param rr vector of scales
       *
       * @param l the order \f$l\f$
       *
       * @param l1 the order \f$l_1\f$
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       */
      void eff_l_l1 (std::vector<std::vector<double>> &eff, const std::vector<double> rr, const int l, const int l1, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the quantity \f$ I_{\mathcal{L} l} (r_1, r_2)\f$
       *
       * This function computes the quantity \f$ I_{\mathcal{L} l} (r_1, r_2)\f$:
       *
       * \f[ I_{\mathcal{L} l} (r_1, r_2) = \sum_{l_1}
       *   (-1)^{l_1+l}(2l_1+1)(2l+1) \begin{pmatrix} l_1 & l &
       *   \mathcal{L} \\ 0 & 0 & 0 \end{pmatrix}^2 \\ \times \int
       *   \mathrm{d} r \;r\;f_{l, l_1}(r_1;r) f_{l, l_1} (r_2; r) \f]
       *
       * where \f$ f_{l, l_1} (r_i;r) \f$ is computed by
       * cbl::cosmology::Cosmology::eff_l_l1 This quantity is used the
       * compute the tree-level theoretical prediction for the biased
       * and redshift space the three-point correlation function,
       * following Slepian&Eisenstein, 2017
       *
       * @param II the quantity \f$ I_{\mathcal{L} l} (r_1, r_2) \f$
       * 
       * @param rr vector of scales
       *
       * @param ll the order \f$l\f$
       *
       * @param LL the order \f$ \mathcal{L} \f$
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       */
      void I_ELL_ell (std::vector<std::vector<double>> &II, const std::vector<double> rr, const int ll, const int LL, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the quantity \f$ k_l (r_1, r_2) \f$
       *
       * This function computes the quantity \f$ k_l (r_1, r_2) \f$:
       *
       * \f[
       *  k_l (r_1, r_2) = \frac{64}{77175} \left[9I_{1,l}(r_1, r_2)-
       *  14I_{3,l}(r_1, r_2) +5I_{5,l}(r_1, r_2) \right] 
       * \f]
       *
       * where \f$ I_{\mathcal{L}l} \f$ is computed by cbl::cosmology::Cosmology::I_ELL_ell
       * This quantity is  used the compute the tree-level theoretical 
       * prediction for the biased and redshift space three-point correlation function,
       * following Slepian&Eisenstein, 2017
       *
       * @param KK the quantity \f$ k_l (r_1, r_2) \f$
       * 
       * @param rr vector of scales
       *
       * @param ll the order \f$l\f$
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       */
      void k_ell (std::vector<std::vector<double>> &KK, const std::vector<double> rr, const int ll, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief the multiplicative factor for \f$ \zeta_0 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_0 \f$, with  local bias:
       *
       * \f[
       *  l = 0 : b_1^4\left[\frac{2}{3}\beta + \frac{38}{45}\beta^2 + \frac{2}{5}\beta^3 + \frac{2}{25}\beta^4\right] + b_1^3\Biggl[\frac{34}{21}\left(1 + \frac{47}{51}\beta + \frac{163}{425}\beta^2 + \frac{201}{2975}\beta^3\right) + \gamma\left(1 + \frac{2}{3}\beta + \frac{1}{9}\beta^2\right) -\frac{4}{3}\gamma_t\left(1 + \frac{2}{3}\beta + \frac{7}{75}\beta^2\right)\Biggr]
       * \f]
       *
       * with \f$b_1\f$ the linear bias, \f$\gamma\f$ the ratio of quadratic and linear bias, 
       * \f$\gamma= 2 b_2 / b_1 \f$ and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param gamma the ratio of quadratic and linear bias, \f$\gamma= 2 b_2 / b_1 \f$, 
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_0 \f$, with local bias
       */
      double zeta_ell_0_factor (const double b1, const double gamma, const double beta) const;


      /**
       * @brief the multiplicative factor for \f$ \zeta_1 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_1 \f$, with local bias:
       *
       * \f[
       *   l = 1 : b_1^4\left[\frac{1}{3}\beta + \frac{3}{5}\beta^2 + \frac{67}{175}\beta^3 + \frac{3}{35}\beta^4\right] + b_1^3\left[1 + \beta + \frac{37}{75}\beta^2 + \frac{17}{175}\beta^3\right]
       * \f]
       *
       * with \f$b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_1 \f$, with local bias
       */
      double zeta_ell_1_factor (const double b1, const double beta) const;

      /**
       * @brief the multiplicative factor for \f$ \zeta_2 \f$, with
       * local bias
       *
       * This function computes the multiplicative factor for \f$
       * \zeta_2 \f$, with local bias:
       *
       * \f[ 
       *     l = 2 : b_1^4\left[\frac{16}{45}\beta^2 + \frac{16}{35}\beta^3 + \frac{32}{245}\beta^4\right] + b_1^3\Biggl[\frac{8}{21}\left(1 + \frac{4}{3}\beta + \frac{54}{35}\beta^2 + \frac{111}{245}\beta^3\right) + \frac{4\beta^2\gamma}{45} + \frac{4}{3}\gamma_t\left(1 + \frac{2}{3}\beta + \frac{1}{21}\beta^2\right)\Biggr]
       * \f]
       *
       * with \f$b_1\f$ the linear bias, \f$\gamma\f$ the ratio of quadratic and linear bias, 
       * \f$\gamma= 2 b_2 / b_1 \f$ and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param gamma the ratio of quadratic and linear bias, \f$\gamma= 2 b_2 / b_1 \f$, 
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_2 \f$, with local bias
       */
      double zeta_ell_2_factor (const double b1, const double gamma, const double beta) const;

      /**
       * @brief the multiplicative factor for \f$ \zeta_3 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_3 \f$, with local bias:
       *
       * \f[
       *    l=3: b_1^4\left[\frac{8}{175}\beta^3 + \frac{8}{315}\beta^4\right] + b_1^3\left[\frac{8}{75}\beta^2 + \frac{8}{175}\beta^3\right]
       * \f]
       *
       * with \f$b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_3 \f$, with local bias
       */
      double zeta_ell_3_factor (const double b1, const double beta) const;

      /**
       * @brief the multiplicative factor for \f$ \zeta_4 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_4 \f$, with local bias:
       *
       * \f[
       *   l=4: b_1^4 \, \frac{128}{11025}\beta^4 + b_1^3 \beta^2\left(-\frac{32}{3675} + \frac{32}{8575}\right)
       * \f]
       *
       * with \f$ b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_4 \f$, with local bias
       */
      double zeta_ell_4_factor (const double b1, const double beta) const;

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l>4 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_k, l>4 \f$, with  local bias:
       *
       * \f[
       *   l>4: b_1^3(7\beta^2+3\beta^3)
       * \f]
       *
       * with \f$b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_l, l>4 \f$, with local bias
       */
      double zeta_ell_k_factor (const double b1, const double beta) const;

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l=0 \f$, with tidal bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_0 \f$, with tidal bias:
       *
       * \f[
       *   l=0: b_1^3\left[\frac{4\gamma_t}{225}\left(-75-50\beta -7\beta^2\right) \right]
       * \f]
       *
       * with \f$\gamma_t\ = b_t/b_1\f$ the ratio between the tidal and the linear bias 
       * and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param gamma_t the ratio between the tidal and the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_l, l=0 \f$, with tidal bias
       */
      double zeta_ell_0_factor_tidal (const double b1, const double gamma_t, const double beta) const;

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l=2 \f$, with tidal bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_2 \f$, with tidal bias:
       *
       * \f[
       *   l=2: b_1^3\left[\frac{4\gamma_t}{63} \left(21 + 14\beta + \beta^2\right)\right]
       * \f]
       *
       * with \f$\gamma_t\ = b_t/b_1\f$ the ratio between the tidal and the linear bias 
       * and \f$ \beta = f/b_1\f$ with \f$ f \f$ the linear growth rate
       *
       * @param gamma_t the ratio between the tidal and the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_l, l=2 \f$, with tidal bias
       */
      double zeta_ell_2_factor_tidal (const double b1, const double gamma_t, const double beta) const;

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l=4 \f$, with tidal bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_4 \f$, with tidal bias:
       *
       * \f[
       *   l=4: b_1^3\,\left[\frac{32 \beta^2 \gamma_t} {525} \right]
       * \f]
       *
       * with \f$\gamma_t\ = b_t/b_1\f$ the ratio between the tidal and the linear bias 
       * and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param gamma_t the ratio between the tidal and the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1\f$ with \f$ f \f$ the linear growth rate
       *
       * @return the multiplicative factor for \f$ \zeta_l, l=4 \f$, with tidal bias
       */
      double zeta_ell_4_factor_tidal (const double b1, const double gamma_t, const double beta) const;

      /**
       *  @brief the pre-cyclic redshift space three-point correlation function as
       *  described in Slepian & Eisenstein 2017
       *
       *  this function computes the redshift space pre-cyclic three-point
       *  correlation function as described in Slepian & Eisenstein 2017, as
       *  follows:
       *
       *  \f[ \zeta_{pc} = \sum_{l=0}^4 \zeta_{pc l}(r_1, r_2)
       *		   P_l(\hat{r_1}\cdot \hat{r_2})+ \sum_{l=0}^4
       *		   \zeta_{pc l}(r_2, r_3) P_l(\hat{r_2}\cdot
       *		   \hat{r_3})+ \sum_{l=0}^4 \zeta_{pc l}(r_3,
       *		   r_1) P_l(\hat{r_3}\cdot \hat{r_1}) \f]
       *
       *  with 
       *
       *  \f[r_3 = \sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos(\theta)}\f]
       *
       *  and
       *
       *  \f[
       *                \zeta_{pc, 0} (\mathbf{r_i}, \mathbf{r_j}) = \quad &\xi^{[0]}(\mathbf{r_i})\xi^{[0]}(\mathbf{r_j})\Biggl\{b_1^4\left[\frac{2}{3}\beta + \frac{38}{45}\beta^2 + \frac{2}{5}\beta^3 + \frac{2}{25}\beta^4\right]\\
    &+ b_1^3\Biggl[\frac{34}{21}\left(1 + \frac{47}{51}\beta + \frac{163}{425}\beta^2 + \frac{201}{2975}\beta^3\right)\\
    &+ \gamma\left(1 + \frac{2}{3}\beta + \frac{1}{9}\beta^2\right) -\frac{4}{3}\gamma_t\left(1 + \frac{2}{3}\beta + \frac{7}{75}\beta^2\right)\Biggr]
    \Biggr\}\\
    &+ b_1^3\beta^2\left(7 + 3\beta\right) \kappa_{0}(\mathbf{r_i}, \mathbf{r_j}),
       *  \f]
       *
       *  \f[
       *                \zeta_{pc, 1} (\mathbf{r_i}, \mathbf{r_j}) = \quad &-\bigg[\xi^{[1+]}(\mathbf{r_i})\xi^{[1-]}(\mathbf{r_j})+\xi^{[1+]}(\mathbf{r_j})\xi^{[1-]}(\mathbf{r_i})\bigg]\Biggl\{b_1^4\left[\frac{1}{3}\beta + \frac{3}{5}\beta^2 + \frac{67}{175}\beta^3 + \frac{3}{35}\beta^4\right]\\
    &+ b_1^3\left[1 + \beta + \frac{37}{75}\beta^2 + \frac{17}{175}\beta^3\right]
    \Biggr\} + b_1^3\beta^2\left(7 + 3\beta\right) \kappa_{1}(\mathbf{r_i}, \mathbf{r_j}),
       *  \f]
       *
       *
       *  \f[
       *                \zeta_{pc, 2} (\mathbf{r_i}, \mathbf{r_j}) = \quad &\xi^{[2]}(\mathbf{r_i})\xi^{[2]}(\mathbf{r_j})\Biggl\{b_1^4\left[\frac{16}{45}\beta^2 + \frac{16}{35}\beta^3 + \frac{32}{245}\beta^4\right]\\
    &+ b_1^3\Biggl[\frac{8}{21}\left(1 + \frac{4}{3}\beta + \frac{54}{35}\beta^2 + \frac{111}{245}\beta^3\right)\\
    &+ \frac{4\beta^2\gamma}{45} + \frac{4}{3}\gamma_t\left(1 + \frac{2}{3}\beta + \frac{1}{21}\beta^2\right)\Biggr]
    \Biggr\} + b_1^3\beta^2\left(7 + 3\beta\right) \kappa_{2}(\mathbf{r_i}, \mathbf{r_j}),
       *  \f]
       *
       *  \f[
       *                \zeta_{pc, 3} (\mathbf{r_i}, \mathbf{r_j}) = \quad &-\bigg[\xi^{[3+]}(\mathbf{r_i})\xi^{[3-]}(\mathbf{r_j})+\xi^{[3+]}(\mathbf{r_j})\xi^{[3-]}(\mathbf{r_i})\bigg]\Biggl\{b_1^4\left[\frac{8}{175}\beta^3 + \frac{8}{315}\beta^4\right]\\
    &+ b_1^3\left[\frac{8}{75}\beta^2 + \frac{8}{175}\beta^3\right]
    \Biggr\} + b_1^3\beta^2\left(7 + 3\beta\right) \kappa_{3}(\mathbf{r_i}, \mathbf{r_j}),
       *  \f]
       *
       *  \f[
       *                \zeta_{pc, 4} (\mathbf{r_i}, \mathbf{r_j}) = \quad &\xi^{[4]}(\mathbf{r_i})\xi^{[4]}(\mathbf{r_j})\Biggl\{\frac{128}{11025}\beta^4b_1^4 + b_1^3\left[-\frac{32}{3675}\beta^2 + \frac{32}{8575}\beta^3 + \frac{32}{525}\beta^{2} \gamma_{t}\right]
    \Biggr\}\\
    &+ b_1^3\beta^2\left(7 + 3\beta\right) \kappa_{4}(\mathbf{r_i}, \mathbf{r_j}),
       *  \f]
       *
       *  where \f$ b_1, b_2, b_t, \beta \f$ are the linear bias, the non-linear bias, the tidal tensor bias and the kaiser factor
       *  respectively, and \f$\xi_{DM}^{[n]}(r)\f$, \f$\xi^{[n\pm]}_{DM}(r)\f$,
       *  are the integrals of the dark matter
       *  power spectrum computed by
       *  cbl::cosmology::ThreePointCorrelation::xi_r_n and cbl::cosmology::ThreePointCorrelation::xi_r_pm.
       *  The coefficients \f$ \zeta_{pc, l} \f$ are computed by:
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_0_factor,
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_1_factor,
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_2_factor,
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_3_factor,
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_4_factor,
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_k_factor,
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_0_factor_tidal, 
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_2_factor_tidal, 
       *  cbl::cosmology::ThreePointCorrelation::zeta_ell_4_factor_tidal.
       *
       *  @param r1 the first side
       *
       *  @param r2 the second side
       *
       *  @param mu the cosine of the angle between r1 and r2
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param bt the tidal tensor bias
       *
       *  @param beta the kaiser factor
       *
       *  @param interp_xil vector containing \f$\xi^{[n]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param interp_xil_m vector containing
       *  \f$\xi^{[n-]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param interp_xil_p vector containing
       *  \f$\xi^{[n+]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param interp_Kl vector containing
       *  \f$\kappa_{\ell}(r_1, r_2)\f$ as defined in Slepian & Eisenstein 2017 (eq. 17)
       *
       *  @param use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return the pre-cyclic redshift space three-point
       *  correlation function
       */
      double zeta_precyclic_RSD (const double r1, const double r2, const double mu, const double b1, const double b2, const double bt, const double beta, std::vector<std::shared_ptr<cbl::glob::FuncGrid>> interp_xil, std::vector<std::shared_ptr<cbl::glob::FuncGrid>> interp_xil_m, std::vector<std::shared_ptr<cbl::glob::FuncGrid>> interp_xil_p, std::vector<std::shared_ptr<cbl::glob::FuncGrid2D>> interp_Kl, const bool use_k) const;

      /**
       *  @brief the terms of the \f$\zeta(r_1, r_2)\f$ expansion 
       *
       *  this function computes the terms of the \f$\zeta(r_1,
       *  r_2)\f$ expansion up to an arbitrary order \f$l\f$ (the
       *  default value is \f$l_{max}=9\f$), as described in Slepian
       *  et al. 2015:
       * 
       *  \f[\zeta_l(r_1, r_2) = \frac{2l+1}{2} \int_{-1}^{1}
       *  \mathrm{d}\mu_{12} \left[\zeta_{pc}(r_1, r_2,
       *  \mu_{12})+\zeta_{pc}(r_2, r_3, \mu_{23})+ \zeta_{pc}(r_3,
       *  r_1, \mu_{31})\right] P_l(\mu_{12}) .\f] 
       * 
       *  the terms in square brackets is computed by
       *  cbl::cosmology::ThreePointCorrelation::zeta_precyclic_RSD
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] b1 the linear bias of the triangle
       *
       *  @param [in] b2 the non-linear bias
       *
       *  @param [in] bt the tidal tensor bias
       *
       *  @param [in] beta the kaiser factor
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xil vector containing \f$\xi^{[n]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param [out] xil_m vector containing
       *  \f$\xi^{[n-]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param [out] xil_p vector containing
       *  \f$\xi^{[n+]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param [out] Kl vector containing
       *  \f$\kappa_{\ell}(r_1, r_2)\f$ as defined in Slepian & Eisenstein 2017 (eq. 17)
       *
       *  @param [in] norders the maximum numbers of orders
       *
       *  @param [in] prec the integral precision
       *
       *  @param [in] use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return vector containing the terms of legendre expansion
       */
      std::vector<double> zeta_expansion_RSD (const double r1, const double r2, const double b1, const double b2, const double bt, const double beta, std::vector<double> &rr, std::vector<std::vector<double>> &xil, std::vector<std::vector<double>> &xil_m, std::vector<std::vector<double>> &xil_p, std::vector<std::vector<std::vector<double>>> &Kl, const int norders=9, const double prec=1.e-3, const bool use_k=false) const;

      /**
       *  @brief the redshift space dark matter three-point correlation function
       *  model by Slepian & Eisenstein 2017
       *
       *  this function computes \f$\zeta_{DM} (r_1, r_2, \hat{r_1}
       *  \cdot \hat{r_2})\f$, as described in Slepian & Eisenstein 2017:
       *
       *  \f[ \zeta_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = 
       *  \sum_l \zeta_l(r_1, r_2) P_l(\hat{r_1} \cdot \hat{r_2}) .\f]
       *
       *  The coefficients of the expansion are computed by
       *  cbl::cosmology::ThreePointCorrelation::zeta_expansion_RSD (using \f$b_1 = 1\f$, \f$b_2 = 0\f$ and \f$b_t = 0\f$)
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       *
       *  @param [in] beta the kaiser factor
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xil vector containing \f$\xi^{[n]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param [out] xil_m vector containing
       *  \f$\xi^{[n-]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param [out] xil_p vector containing
       *  \f$\xi^{[n+]}_{DM}(r)\f$ as defined in Slepian & Eisenstein 2017 (eq. 20)
       *
       *  @param [out] Kl vector containing
       *  \f$\kappa_{\ell}(r_1, r_2)\f$ as defined in Slepian & Eisenstein 2017 (eq. 17)
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_matter the dark matter power spectrum
       *
       *  @param [in] norders the maximum number of orders
       *
       *  @param [in] prec the integral precision
       *
       *  @param [in] use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return the redshift space connected dark matter three-point correlation
       *  function
       */
      double zeta_DM_Slepian_RSD (const double r1, const double r2, const double theta, const double beta, std::vector<double> &rr, std::vector<std::vector<double>> &xil, std::vector<std::vector<double>> &xil_m, std::vector<std::vector<double>> &xil_p, std::vector<std::vector<std::vector<double>>> &Kl, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders=9, const double prec=1.e-3, const bool use_k=false);


      /**
       *  @brief the local-bias model of the three-point correlation
       *  function of dark matter haloes in redshift space
       *
       *  this function computes the three-point correlation function
       *  of dark matter haloes in redshift space with the Slepian & Eisenstein (2017)
       *  model, as follows:
       *
       *  \f[ \zeta_h (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = b_1^3
       *  \zeta_{DM}(r_1, r_2, \hat{r_1} \cdot \hat{r_2}) + b_1^2 b_2
       *  \left[ \xi(r_1)\cdot\xi(r_2) + \xi(r_2)\cdot\xi(r_3) +
       *  \xi(r_3)\cdot\xi(r_1) \right] \f]
       *
       *  with \f$r_3 = \sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos(\theta)}\f$
       *  and \f$b_1, b_2\f$ the linear and non-linear halo bias,
       *  respectively; \f$\zeta_{DM}\f$ is computed by
       *  cbl::cosmology::ThreePointCorrelation::zeta_DM_Slepian_RSD
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *  
       *  @param theta the vector of angles between r1 and r2
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param bt the tidal tensor bias
       *
       *  @param beta the kaiser factor
       *
       *  @param model the model to compute the three-point
       *  correlation function. For now, it can only be "Slepian"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @param norders the maximum number of orders
       *
       *  @param use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return vector containing the three-point correlation
       *  function of dark matter haloes in redshift space
       */
      std::vector<double> zeta_halo_RSD (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double bt, const double beta, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders = 9, const bool use_k=false);


      /**
       *  @brief the non-local-bias model of the three-point correlation
       *  function of galaxies in redshift space
       *
       *  this function computes the three-point correlation function
       *  of galaxies in redshift space with the Slepian & Eisenstein (2017)
       *  model, as follows:
       *
       *  \f[ \zeta_{g} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = 
       *  \sum_l \zeta_l(r_1, r_2) P_l(\hat{r_1} \cdot \hat{r_2}) .\f]
       *
       *  The coefficients \f[ \zeta_l(r_1, r_2) \f] are computed using the function cbl::cosmology::ThreePointCorrelation::zeta_expansion_RSD
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *  
       *  @param theta the angle between r1 and r2
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param bt the tidal tensor bias
       *
       *  @param beta the kaiser factor
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @param norders the maximum number of orders
       *
       *  @param prec the integral precision
       *
       *  @param use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return three-point correlation
       *  function of galaxies in redshift space
       */
      double zeta_RSD (const double r1, const double r2, const double theta, const double b1, const double b2, const double bt, const double beta, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders=9, const double prec=1.e-3, const bool use_k=false);

      /**
       *  @brief the non-local-bias model of the three-point correlation
       *  function of galaxies in redshift space
       *
       *  this function computes the three-point correlation function
       *  of galaxies in redshift space with the Slepian & Eisenstein (2017)
       *  model, as follows:
       *
       *  \f[ \zeta_{g} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = 
       *  \sum_l \zeta_l(r_1, r_2) P_l(\hat{r_1} \cdot \hat{r_2}) .\f]
       *
       *  The coefficients \f[ \zeta_l(r_1, r_2) \f] are computed using the function cbl::cosmology::ThreePointCorrelation::zeta_expansion_RSD
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *  
       *  @param theta the vector of angles between r1 and r2
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param bt the tidal tensor bias
       *
       *  @param beta the kaiser factor
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @param norders the maximum number of orders
       *
       *  @param prec the integral precision
       *
       *  @param use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return vector containing the three-point correlation
       *  function of galaxies in redshift space
       */
      std::vector<double> zeta_RSD (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double bt, const double beta, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders=9, const double prec=1.e-3, const bool use_k=false);

      /**
       *  @brief the non-local-bias model of the three-point correlation
       *  function of galaxies in redshift space at all scales
       *
       *  this function computes the three-point correlation function
       *  of galaxies in redshift space at all scales. The three-point correlation
       *  function at fixed scales is computed using the function cbl::cosmology::ThreePointCorrelation::zeta_RSD.
       *
       *  @param index vector containing triangle indices
       *
       *  @param r12 vector containing the first sides of the triangles
       *
       *  @param r13 vector containing the second sides of the triangles
       *  
       *  @param r23 vector containing the third sides of the triangles
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param bt the tidal tensor bias
       *
       *  @param beta the kaiser factor
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_matter the dark matter power spectrum
       *
       *  @param norders the maximum number of orders
       *
       *  @param prec the integral precision
       *
       *  @param use_k if true, use the \f$\kappa_{\ell}\f$ part of the model
       * \f$O(\beta^2)\f$
       *
       *  @return vector containing the three-point correlation
       *  function of galaxies in redshift space at all scales
       */
      std::vector<double> zeta_RSD_all_scales (const std::vector<double> index, const std::vector<double> r12, const std::vector<double> r13, const std::vector<double> r23, const double b1, const double b2, const double bt, const double beta, const std::vector<double> kk, const std::vector<double> Pk_matter, const int norders=9, const double prec=1.e-3, const bool use_k=false);

      ///@}
      
    };
  }
}

#endif
