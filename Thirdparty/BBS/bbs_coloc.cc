/**
* This file is part of DefSLAM.
*
* Copyright (C) 2017-2020 Jose Lamarca Peiro <jlamarca at unizar dot es>, J.M.M. Montiel (University
*of Zaragoza) && Shaifali Parashar, Adrien Bartoli (Université Clermont Auvergne)
* For more information see <https://github.com/unizar/DefSLAM>
*
* DefSLAM is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DefSLAM is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DefSLAM. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NTHREADS
#define NTHREADS 2
#endif
#include <bbs_coloc.h>
#include <iostream>
#include <bbs_MAC.h>

#define max(a, b) (((a) >= (b)) ? (a) : (b))
#define min(a, b) (((a) <= (b)) ? (a) : (b))
namespace BBS
{
  /* Build a sparse colocation matrix.
 * This function does not handle sites that are outside of the domain (an error
 * is issued in such cases).
 * ARGUMENTS:
 * - bbs: b-spline structure.
 * - u,v: sites used to build the colocation matrix.
 * - nsites: number of sites (number of elements in u and v).
 * - pr: array containing the values of the matrix (size: number of non-zero
 * elements, i.e. 16*nsites).
 * - ir: ir(i) is the row index (0-based) of pr(i). The sizes of ir and pr are
 * equal.
 * - jc: jc(j) is the shift in pr (and ir) for the j-th column. numel(jc) =
 * 1+nptsx*nptsy (number of columns).
 * RETURN VALUE:
 * This function return 0 if successful, an error code otherwise:
 * 0: successful.
 * 1: a colocation site was ouside of the spline definition domain. */
  typedef void (*basis_func_t)(double, double *);

  /* Return a pointer to the basis function that correspond to the "order" of
 * derivation. */
  basis_func_t get_basis_ptr_e(int order)
  {
    switch (order)
    {
    case 0:
      return eval_basis;
    case 1:
      return eval_basis_d;
    case 2:
      return eval_basis_dd;
    }
    return nullptr;
  }

  double get_deriv_fact_e(bbs_t *bbs, int uorder, int vorder)
  {
    double su = (bbs->umax - bbs->umin) / (bbs->nptsu - 3);
    double sv = (bbs->vmax - bbs->vmin) / (bbs->nptsv - 3);

    return 1.0 / (pow(su, uorder) * pow(sv, vorder));
  }

  void colocEigen(bbs_t *bbs, double *u, double *v, int nsites,
                  Eigen::MatrixXd &colocMatrix)
  {
    int k, iu, iv, col, Iu, Iv;
    double nu[nsites];
    double nv[nsites];
    int interu[nsites];
    int interv[nsites];
    double basis_u[4];
    double basis_v[4];
    int npts = bbs->nptsu * bbs->nptsv;

    size_t nb_elem_col[npts];
    size_t cur_ind[npts];
    size_t jc[npts];
    jc[0] = 0;
    // Compute the normalized evaluation values and their interval numbers
    std::cout << bbs->umin << "  " << bbs->umax << " " << bbs->vmin;

    normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
    normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);

    // FILLING jc
    // First, compute the number of elements in each column
    // Check at the same time if the sites are inside the definition domain.
    for (k = 0; k < nsites; ++k)
    {
      Iu = interu[k];
      Iv = interv[k];

      if (Iu < 0 || Iu > bbs->nptsu - 4 || Iv < 0 || Iv > bbs->nptsv - 4)
      {
        return;
      }

      for (iu = 0; iu < 4; ++iu)
        for (iv = 0; iv < 4; ++iv)
          nb_elem_col[(iu + Iu) * bbs->nptsv + iv + Iv] += 1;
    }

    // Second, compute jc from nb_elem_col
    for (k = 1; k <= npts; ++k)
      jc[k] = jc[k - 1] + nb_elem_col[k - 1];

    // FILLING pr and ir
    // cur_ind contains the current shift in pr (and ir) for each column
    for (k = 0; k < npts; ++k)
      cur_ind[k] = jc[k];

    for (k = 0; k < nsites; ++k)
    {
      eval_basis(nu[k], basis_u);
      eval_basis(nv[k], basis_v);
      for (iu = 0; iu < 4; ++iu)
        for (iv = 0; iv < 4; ++iv)
        {
          col = (iu + interu[k]) * bbs->nptsv + iv + interv[k];
          colocMatrix(k, col) = basis_u[iu] * basis_v[iv];
          cur_ind[col] += 1;
        }
    }
  }

  int coloc_derivEigen(bbs_t *bbs, double *u, double *v, int nsites, int du,
                       int dv, Eigen::MatrixXd &colocMatrix)
  {
    int k, iu, iv, col, Iu, Iv;
    int ret_code = 0;
    double nu[nsites];
    double nv[nsites];
    int interu[nsites];
    int interv[nsites];
    basis_func_t b_func_u = get_basis_ptr_e(du);
    basis_func_t b_func_v = get_basis_ptr_e(dv);
    double fact = get_deriv_fact_e(bbs, du, dv);
    double basis_u[4];
    double basis_v[4];
    int npts = bbs->nptsu * bbs->nptsv;
    size_t nb_elem_col[npts];
    size_t cur_ind[npts];
    double jc[npts];
    jc[0] = 0;
    nb_elem_col[0] = 0;
    // double ir[npts];
    // Compute the normalized evaluation values and their interval numbers
    normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
    normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);
    //
    // FILLING jc
    // First, compute the number of elements in each column
    // Check at the same time if the sites are inside the definition domain.
    for (k = 0; k < nsites; ++k)
    {
      Iu = interu[k];
      Iv = interv[k];

      if (Iu < 0 || Iu > bbs->nptsu - 4 || Iv < 0 || Iv > bbs->nptsv - 4)
      {
        ret_code = 1;
        goto error;
      }

      for (iu = 0; iu < 4; ++iu)
        for (iv = 0; iv < 4; ++iv)
          nb_elem_col[(iu + Iu) * bbs->nptsv + iv + Iv] += 1;
    }

    // Second, compute jc from nb_elem_col
    for (k = 1; k <= npts; ++k)
      jc[k] = jc[k - 1] + nb_elem_col[k - 1];

    // FILLING pr and ir
    // cur_ind contains the current shift in pr (and ir) for each column
    for (k = 0; k < npts; ++k)
      cur_ind[k] = jc[k];

    for (k = 0; k < nsites; ++k)
    {
      b_func_u(nu[k], basis_u);
      b_func_v(nv[k], basis_v);
      for (iu = 0; iu < 4; ++iu)
        for (iv = 0; iv < 4; ++iv)
        {
          col = (iu + interu[k]) * bbs->nptsv + iv + interv[k];
          colocMatrix(k, col) = fact * basis_u[iu] * basis_v[iv];
          cur_ind[col] += 1;
        }
    }

  error:
    return ret_code;
  }

  using namespace std;
  double __bxx_coeff1[] = {
      1.0 / 756.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 504.0, 1.0 / 252.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, -1.0 / 504.0, 1.0 / 252.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 1512.0, 0.0, -1.0 / 504.0, 1.0 / 756.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      43.0 / 5040.0, -43.0 / 3360.0, 0.0, 43.0 / 10080.0,
      11.0 / 140.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0, 0.0,
      -33.0 / 280.0, 33.0 / 140.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0,
      0.0, -33.0 / 280.0, 33.0 / 140.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      43.0 / 10080.0, 0.0, -43.0 / 3360.0, 43.0 / 5040.0,
      11.0 / 280.0, 0.0, -33.0 / 280.0, 11.0 / 140.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 252.0, -1.0 / 168.0, 0.0, 1.0 / 504.0,
      311.0 / 5040.0, -311.0 / 3360.0, 0.0, 311.0 / 10080.0,
      11.0 / 140.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0, 0.0,
      -311.0 / 3360.0, 311.0 / 1680.0, -311.0 / 3360.0, 0.0,
      -33.0 / 280.0, 33.0 / 140.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0,
      0.0, -311.0 / 3360.0, 311.0 / 1680.0, -311.0 / 3360.0,
      0.0, -33.0 / 280.0, 33.0 / 140.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 504.0, 0.0, -1.0 / 168.0, 1.0 / 252.0,
      311.0 / 10080.0, 0.0, -311.0 / 3360.0, 311.0 / 5040.0,
      11.0 / 280.0, 0.0, -33.0 / 280.0, 11.0 / 140.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 15120.0, -1.0 / 10080.0, 0.0, 1.0 / 30240.0,
      1.0 / 252.0, -1.0 / 168.0, 0.0, 1.0 / 504.0,
      43.0 / 5040.0, -43.0 / 3360.0, 0.0, 43.0 / 10080.0,
      1.0 / 756.0, 0.0, 0.0, 0.0,
      -1.0 / 10080.0, 1.0 / 5040.0, -1.0 / 10080.0, 0.0,
      -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0, 0.0,
      -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0, 0.0,
      -1.0 / 504.0, 1.0 / 252.0, 0.0, 0.0,
      0.0, -1.0 / 10080.0, 1.0 / 5040.0, -1.0 / 10080.0,
      0.0, -1.0 / 168.0, 1.0 / 84.0, -1.0 / 168.0,
      0.0, -43.0 / 3360.0, 43.0 / 1680.0, -43.0 / 3360.0,
      0.0, -1.0 / 504.0, 1.0 / 252.0, 0.0,
      1.0 / 30240.0, 0.0, -1.0 / 10080.0, 1.0 / 15120.0,
      1.0 / 504.0, 0.0, -1.0 / 168.0, 1.0 / 252.0,
      43.0 / 10080.0, 0.0, -43.0 / 3360.0, 43.0 / 5040.0,
      1.0 / 1512.0, 0.0, -1.0 / 504.0, 1.0 / 756.0};
  double __bxy_coeff1[] = {
      1.0 / 200.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      7.0 / 1200.0, 17.0 / 600.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 100.0, -29.0 / 1200.0, 17.0 / 600.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 1200.0, -1.0 / 100.0, 7.0 / 1200.0, 1.0 / 200.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      7.0 / 1200.0, 49.0 / 7200.0, -7.0 / 600.0, -7.0 / 7200.0,
      17.0 / 600.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      49.0 / 7200.0, 119.0 / 3600.0, -203.0 / 7200.0, -7.0 / 600.0,
      119.0 / 3600.0, 289.0 / 1800.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -7.0 / 600.0, -203.0 / 7200.0, 119.0 / 3600.0, 49.0 / 7200.0,
      -17.0 / 300.0, -493.0 / 3600.0, 289.0 / 1800.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -7.0 / 7200.0, -7.0 / 600.0, 49.0 / 7200.0, 7.0 / 1200.0,
      -17.0 / 3600.0, -17.0 / 300.0, 119.0 / 3600.0, 17.0 / 600.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 100.0, -7.0 / 600.0, 1.0 / 50.0, 1.0 / 600.0,
      -29.0 / 1200.0, -203.0 / 7200.0, 29.0 / 600.0, 29.0 / 7200.0,
      17.0 / 600.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -7.0 / 600.0, -17.0 / 300.0, 29.0 / 600.0, 1.0 / 50.0,
      -203.0 / 7200.0, -493.0 / 3600.0, 841.0 / 7200.0, 29.0 / 600.0,
      119.0 / 3600.0, 289.0 / 1800.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 50.0, 29.0 / 600.0, -17.0 / 300.0, -7.0 / 600.0,
      29.0 / 600.0, 841.0 / 7200.0, -493.0 / 3600.0, -203.0 / 7200.0,
      -17.0 / 300.0, -493.0 / 3600.0, 289.0 / 1800.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 600.0, 1.0 / 50.0, -7.0 / 600.0, -1.0 / 100.0,
      29.0 / 7200.0, 29.0 / 600.0, -203.0 / 7200.0, -29.0 / 1200.0,
      -17.0 / 3600.0, -17.0 / 300.0, 119.0 / 3600.0, 17.0 / 600.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 1200.0, -7.0 / 7200.0, 1.0 / 600.0, 1.0 / 7200.0,
      -1.0 / 100.0, -7.0 / 600.0, 1.0 / 50.0, 1.0 / 600.0,
      7.0 / 1200.0, 49.0 / 7200.0, -7.0 / 600.0, -7.0 / 7200.0,
      1.0 / 200.0, 0.0, 0.0, 0.0,
      -7.0 / 7200.0, -17.0 / 3600.0, 29.0 / 7200.0, 1.0 / 600.0,
      -7.0 / 600.0, -17.0 / 300.0, 29.0 / 600.0, 1.0 / 50.0,
      49.0 / 7200.0, 119.0 / 3600.0, -203.0 / 7200.0, -7.0 / 600.0,
      7.0 / 1200.0, 17.0 / 600.0, 0.0, 0.0,
      1.0 / 600.0, 29.0 / 7200.0, -17.0 / 3600.0, -7.0 / 7200.0,
      1.0 / 50.0, 29.0 / 600.0, -17.0 / 300.0, -7.0 / 600.0,
      -7.0 / 600.0, -203.0 / 7200.0, 119.0 / 3600.0, 49.0 / 7200.0,
      -1.0 / 100.0, -29.0 / 1200.0, 17.0 / 600.0, 0.0,
      1.0 / 7200.0, 1.0 / 600.0, -7.0 / 7200.0, -1.0 / 1200.0,
      1.0 / 600.0, 1.0 / 50.0, -7.0 / 600.0, -1.0 / 100.0,
      -7.0 / 7200.0, -7.0 / 600.0, 49.0 / 7200.0, 7.0 / 1200.0,
      -1.0 / 1200.0, -1.0 / 100.0, 7.0 / 1200.0, 1.0 / 200.0};
  double __byy_coeff1[] = {
      1.0 / 756.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      43.0 / 5040.0, 11.0 / 140.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 252.0, 311.0 / 5040.0, 11.0 / 140.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 15120.0, 1.0 / 252.0, 43.0 / 5040.0, 1.0 / 756.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 504.0, -43.0 / 3360.0, -1.0 / 168.0, -1.0 / 10080.0,
      1.0 / 252.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -43.0 / 3360.0, -33.0 / 280.0, -311.0 / 3360.0, -1.0 / 168.0,
      43.0 / 1680.0, 33.0 / 140.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 168.0, -311.0 / 3360.0, -33.0 / 280.0, -43.0 / 3360.0,
      1.0 / 84.0, 311.0 / 1680.0, 33.0 / 140.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 10080.0, -1.0 / 168.0, -43.0 / 3360.0, -1.0 / 504.0,
      1.0 / 5040.0, 1.0 / 84.0, 43.0 / 1680.0, 1.0 / 252.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 504.0, -43.0 / 3360.0, -1.0 / 168.0, -1.0 / 10080.0,
      1.0 / 252.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -43.0 / 3360.0, -33.0 / 280.0, -311.0 / 3360.0, -1.0 / 168.0,
      43.0 / 1680.0, 33.0 / 140.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 168.0, -311.0 / 3360.0, -33.0 / 280.0, -43.0 / 3360.0,
      1.0 / 84.0, 311.0 / 1680.0, 33.0 / 140.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 10080.0, -1.0 / 168.0, -43.0 / 3360.0, -1.0 / 504.0,
      1.0 / 5040.0, 1.0 / 84.0, 43.0 / 1680.0, 1.0 / 252.0,
      0.0, 0.0, 0.0, 0.0,
      1.0 / 1512.0, 43.0 / 10080.0, 1.0 / 504.0, 1.0 / 30240.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 504.0, -43.0 / 3360.0, -1.0 / 168.0, -1.0 / 10080.0,
      1.0 / 756.0, 0.0, 0.0, 0.0,
      43.0 / 10080.0, 11.0 / 280.0, 311.0 / 10080.0, 1.0 / 504.0,
      0.0, 0.0, 0.0, 0.0,
      -43.0 / 3360.0, -33.0 / 280.0, -311.0 / 3360.0, -1.0 / 168.0,
      43.0 / 5040.0, 11.0 / 140.0, 0.0, 0.0,
      1.0 / 504.0, 311.0 / 10080.0, 11.0 / 280.0, 43.0 / 10080.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 168.0, -311.0 / 3360.0, -33.0 / 280.0, -43.0 / 3360.0,
      1.0 / 252.0, 311.0 / 5040.0, 11.0 / 140.0, 0.0,
      1.0 / 30240.0, 1.0 / 504.0, 43.0 / 10080.0, 1.0 / 1512.0,
      0.0, 0.0, 0.0, 0.0,
      -1.0 / 10080.0, -1.0 / 168.0, -43.0 / 3360.0, -1.0 / 504.0,
      1.0 / 15120.0, 1.0 / 252.0, 43.0 / 5040.0, 1.0 / 756.0};

  void BendingEigen(bbs_t *bbs, double err,
                    Eigen::SparseMatrix<double> &benMatrix)
  {
    int ny = bbs->nptsu; // Order of the dimensions are switched due to a previous
                         // version of the code
    int nx = bbs->nptsv;
    int i(0), j(0), a(0), b(0), c(0), d(0), e1(0), f1(0), e2(0), f2(0), curnb(0),
        total(0), ind(0), sb(0), nb(0);
    double pr[bbs->nptsu * bbs->nptsv * bbs->nptsu * bbs->nptsv];
    size_t ir[bbs->nptsu * bbs->nptsv];
    size_t jc[bbs->nptsu * bbs->nptsv];

    size_t *pi = ir, *pp = jc;
    double *px = pr;
    double coeff_b[256];
    double lbd;
    double sy = (bbs->umax - bbs->umin) / (bbs->nptsu - 3);
    double sx = (bbs->vmax - bbs->vmin) / (bbs->nptsv - 3);
    *pp = 0;
    ++pp;
    total = 0;
    /* Location of the non-zeros coefficients of the matrix. */
    for (j = 0; j < ny; ++j)
    {
      for (i = 0; i < nx; ++i)
      {
        curnb = 0;
        for (b = max(j - 3, 0); b <= j - 1; ++b)
          for (a = max(0, i - 3); a <= min(nx - 1, i + 3); ++a)
          {
            *pi = b * nx + a;
            *px = 0.0;
            ++pi;
            ++px;
            ++curnb;
          }

        for (a = max(0, i - 3); a <= i; ++a)
        {
          *pi = j * nx + a;
          *px = 0.0;
          ++pi;
          ++px;
          ++curnb;
        }
        total += curnb;
        *pp = total;
        ++pp;
      }
    }

    // Build the coefficients of the B matrix with the scales taken into account
    for (j = 0; j < 16; ++j)
      for (i = 0; i <= j; ++i)
      {
        ind = 16 * j + i;
        coeff_b[ind] = sy * __bxx_coeff1[ind] / pow(sx, 3) +
                       __bxy_coeff1[ind] / (sx * sy) +
                       sx * __byy_coeff1[ind] / pow(sy, 3);
      }

    // Put the right coefficients at the right locations (one knot domain at a
    // time
    // and with the scaling given by lambda)
    pp = jc;
    px = pr;
    for (b = 0; b < ny - 3; ++b)
    {
      for (a = 0; a < nx - 3; ++a)
      {
        lbd = err;
        for (c = 0; c < 16; ++c)
        {
          for (d = c; d < 16; ++d)
          {
            e1 = c / 4;
            f1 = c % 4;
            e2 = d / 4;
            f2 = d % 4;
            i = (b + e1) * nx + a + f1;
            j = (b + e2) * nx + a + f2;
            nb = i / nx - max(j / nx - 3, 0);
            sb = min(min(4 + (j % nx), 3 + nx - (j % nx)), min(nx, 7));
            px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)] +=
                lbd * coeff_b[16 * d + c];
            if (i == j)
            {
              benMatrix.coeffRef(i, i) =
                  px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)];
            }
            else
            {
              benMatrix.coeffRef(i, j) =
                  px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)];
              benMatrix.coeffRef(j, i) =
                  px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)];
            }
          }
        }
      }
    }
  }

  void BendingEigen(bbs_t *bbs, double err, Eigen::MatrixXd &benMatrix)
  {
    int ny = bbs->nptsu; // Order of the dimensions are switched due to a previous
                         // version of the code
    int nx = bbs->nptsv;
    int i, j, a, b, c, d, e1, f1, e2, f2, curnb, total, ind, sb, nb;
    double pr[bbs->nptsu * bbs->nptsv * bbs->nptsu * bbs->nptsv];
    size_t ir[bbs->nptsu * bbs->nptsv];
    size_t jc[bbs->nptsu * bbs->nptsv];

    size_t *pi = ir, *pp = jc;
    double *px = pr;
    double coeff_b[256];
    double lbd;
    double sy = (bbs->umax - bbs->umin) / (bbs->nptsu - 3);
    double sx = (bbs->vmax - bbs->vmin) / (bbs->nptsv - 3);
    *pp = 0;
    ++pp;
    total = 0;
    /* Location of the non-zeros coefficients of the matrix. */
    for (j = 0; j < ny; ++j)
    {
      for (i = 0; i < nx; ++i)
      {
        curnb = 0;
        for (b = max(j - 3, 0); b <= j - 1; ++b)
          for (a = max(0, i - 3); a <= min(nx - 1, i + 3); ++a)
          {
            *pi = b * nx + a;
            *px = 0.0;
            ++pi;
            ++px;
            ++curnb;
          }

        for (a = max(0, i - 3); a <= i; ++a)
        {
          *pi = j * nx + a;
          *px = 0.0;
          ++pi;
          ++px;
          ++curnb;
        }
        total += curnb;
        *pp = total;
        ++pp;
      }
    }

    // Build the coefficients of the B matrix with the scales taken into account
    for (j = 0; j < 16; ++j)
      for (i = 0; i <= j; ++i)
      {
        ind = 16 * j + i;
        coeff_b[ind] = sy * __bxx_coeff1[ind] / pow(sx, 3) +
                       __bxy_coeff1[ind] / (sx * sy) +
                       sx * __byy_coeff1[ind] / pow(sy, 3);
      }

    // Put the right coefficients at the right locations (one knot domain at a
    // time
    // and with the scaling given by lambda)
    pp = jc;
    px = pr;
    for (b = 0; b < ny - 3; ++b)
    {
      for (a = 0; a < nx - 3; ++a)
      {
        lbd = err;
        for (c = 0; c < 16; ++c)
        {
          for (d = c; d < 16; ++d)
          {
            e1 = c / 4;
            f1 = c % 4;
            e2 = d / 4;
            f2 = d % 4;
            i = (b + e1) * nx + a + f1;
            j = (b + e2) * nx + a + f2;
            nb = i / nx - max(j / nx - 3, 0);
            sb = min(min(4 + (j % nx), 3 + nx - (j % nx)), min(nx, 7));
            px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)] +=
                lbd * coeff_b[16 * d + c];
            if (i == j)
            {
              benMatrix.coeffRef(i, i) =
                  px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)];
            }
            else
            {
              benMatrix.coeffRef(i, j) =
                  px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)];
              benMatrix.coeffRef(j, i) =
                  px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)];
            }
          }
        }
      }
    }
  }
  /* Pointer type for a basis function evaluation */
  void EvalEigen(bbs_t *bbs, const double *ctrlpts, double *u, double *v,
                 int nsites, Eigen::MatrixXd &val, int du, int dv)
  {
    double nu[nsites];
    double nv[nsites];
    int interu[nsites];
    int interv[nsites];
    // Compute the normalized evaluation values and their interval numbers
    normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
    normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);

#pragma omp parallel for num_threads(NTHREADS), schedule(guided)
    for (int k = 0; k < nsites; ++k)
    {
      int iu, iv, d, ind;
      double basis_u[4];
      double basis_v[4];
      double bu, bas;
      basis_func_t b_func_u = get_basis_ptr_e(du);
      basis_func_t b_func_v = get_basis_ptr_e(dv);
      double fact = get_deriv_fact_e(bbs, du, dv);

      b_func_u(nu[k], basis_u);
      b_func_v(nv[k], basis_v);

      for (d = 0; d < bbs->valdim; ++d)
        val(k, d) = 0.0;

      for (iu = 0; iu < 4; ++iu)
      {
        bu = basis_u[iu];
        for (iv = 0; iv < 4; ++iv)
        {
          bas = bu * basis_v[iv];
          ind = bbs->valdim * ((iu + interu[k]) * bbs->nptsv + iv + interv[k]);
          for (d = 0; d < bbs->valdim; ++d)
            val(k, d) = val(k, d) + ctrlpts[ind++] * bas;
        }
      }

      for (d = 0; d < bbs->valdim; ++d)
        val(k, d) = val(k, d) * fact;
    }
  }
} // namespace BBS
