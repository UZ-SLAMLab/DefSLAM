/* bbs.cpp
 *
 * General routines for cubic bidimensional b-splines.
 *
 * History
 *  2009/??/??: First version
 *  2010/12/15: Translation to CPP
 *              Adding OpenMP support
 *  2011/02/03: Defining number of thread by compilation directive
 *
 * (c)2009-2011, Florent Brunet.
 */

/*
* BBS is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* BBS is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
* for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include "bbs.h"

#ifndef NTHREADS
#define NTHREADS 3
#endif

#define max(a, b) (((a) >= (b)) ? (a) : (b))
#define min(a, b) (((a) <= (b)) ? (a) : (b))
namespace BBS {

/* Normalize the values in "x". "nb_x" is the size of the "x" and "nx" arrays.
 */
void normalize(double xmin, double xmax, int npts, double *x, int nb_x,
               double *nx) {
  int ninter = npts - 3;
  double width_inter = (xmax - xmin) / ninter;

#pragma omp parallel for num_threads(NTHREADS), schedule(guided)
  for (int i = 0; i < nb_x; ++i) {
    if (x[i] == xmax)
      nx[i] = 1.0;
    else if (x[i] < xmin)
      nx[i] = (x[i] - xmin) / width_inter;
    else if (x[i] > xmax)
      nx[i] = (x[i] - xmin) / width_inter - ninter;
    else {
      double scaled = (x[i] - xmin) / width_inter;
      nx[i] = scaled - floor(scaled);
    }
  }
}

/* Same as "normalize" but computes also the interval number to which the x's
 * belongs. */
void normalize_with_inter(double xmin, double xmax, int npts, double *x,
                          int nb_x, double *nx, int *inter) {
  int ninter = npts - 3;
  double width_inter = (xmax - xmin) / ninter;

  //#pragma omp parallel for num_threads(NTHREADS), schedule(guided)
  for (int i = 0; i < nb_x; ++i) {
    if (x[i] == xmax) {
      nx[i] = 1.0;
      inter[i] = ninter - 1;
    } else if (x[i] < xmin) {
      nx[i] = (x[i] - xmin) / width_inter;
      inter[i] = -1;
    } else if (x[i] > xmax) {
      nx[i] = (x[i] - xmin) / width_inter - ninter;
      inter[i] = ninter;
    } else {
      double scaled = (x[i] - xmin) / width_inter;
      inter[i] = (int)floor(scaled);
      nx[i] = scaled - inter[i];
    }
  }
}

/* Evaluate the 4 basis at nx. "nx" must be normalized values. */
void eval_basis(double nx, double *val_basis) {
  double nx2 = nx * nx;
  double nx3 = nx2 * nx;
  val_basis[0] = (-nx3 + 3.0 * nx2 - 3.0 * nx + 1.0) / 6.0;
  val_basis[1] = (3.0 * nx3 - 6.0 * nx2 + 4.0) / 6.0;
  val_basis[2] = (-3.0 * nx3 + 3.0 * nx2 + 3.0 * nx + 1.0) / 6.0;
  val_basis[3] = nx3 / 6.0;
}

/* Evaluate the 4 basis first derivatives at nx. "nx" must be normalized values.
 */
void eval_basis_d(double nx, double *val_basis) {
  double nx2 = nx * nx;
  val_basis[0] = (-nx2 + 2 * nx - 1) / 2.0;
  val_basis[1] = (3.0 * nx2 - 4.0 * nx) / 2.0;
  val_basis[2] = (-3 * nx2 + 2 * nx + 1) / 2.0;
  val_basis[3] = nx2 / 2.0;
}

/* Evaluate the 4 basis second derivatives at nx. "nx" must be normalized
 * values. */
void eval_basis_dd(double nx, double *val_basis) {
  val_basis[0] = -nx + 1.0;
  val_basis[1] = 3.0 * nx - 2.0;
  val_basis[2] = -3.0 * nx + 1.0;
  val_basis[3] = nx;
}

/* Pointer type for a basis function evaluation */
typedef void (*basis_func_t)(double, double *);

/* Return a pointer to the basis function that correspond to the "order" of
 * derivation. */
basis_func_t get_basis_ptr(int order) {
  switch (order) {
  case 0:
    return eval_basis;
  case 1:
    return eval_basis_d;
  case 2:
    return eval_basis_dd;
  }
  return NULL;
}

double get_deriv_fact(bbs_t *bbs, int uorder, int vorder) {
  double su = (bbs->umax - bbs->umin) / (bbs->nptsu - 3);
  double sv = (bbs->vmax - bbs->vmin) / (bbs->nptsv - 3);

  return 1.0 / (pow(su, uorder) * pow(sv, vorder));
}

/* Evaluate a spline (or its derivatives, according to the basis functions
 * utilized).
 * - bbs: bidimensional bspline structure.
 * - ctrlpts: control points (valdim x npts matrix).
 * - u, v: locations where the b-spline is evaluated.
 * - nval: number of values in u (and v).
 * - val: valdim x npts matrix containing the result.
 * - du, dv: order of derivation. */
void eval(bbs_t *bbs, double *ctrlpts, double *u, double *v, int nsites,
          double *val, int du, int dv) {
  double nu[nsites];
  double nv[nsites];
  int interu[nsites];
  int interv[nsites];

  // Compute the normalized evaluation values and their interval numbers
  normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
  normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);

#pragma omp parallel for num_threads(NTHREADS), schedule(guided)
  for (int k = 0; k < nsites; ++k) {
    int iu, iv, d, ind;
    double basis_u[4];
    double basis_v[4];
    double bu, bas;
    basis_func_t b_func_u = get_basis_ptr(du);
    basis_func_t b_func_v = get_basis_ptr(dv);
    double fact = get_deriv_fact(bbs, du, dv);

    b_func_u(nu[k], basis_u);
    b_func_v(nv[k], basis_v);

    for (d = 0; d < bbs->valdim; ++d)
      val[bbs->valdim * k + d] = 0.0;

    for (iu = 0; iu < 4; ++iu) {
      bu = basis_u[iu];
      for (iv = 0; iv < 4; ++iv) {
        bas = bu * basis_v[iv];
        ind = bbs->valdim * ((iu + interu[k]) * bbs->nptsv + iv + interv[k]);
        for (d = 0; d < bbs->valdim; ++d)
          val[bbs->valdim * k + d] += ctrlpts[ind++] * bas;
      }
    }

    for (d = 0; d < bbs->valdim; ++d)
      val[bbs->valdim * k + d] *= fact;
  }
}

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
int coloc(bbs_t *bbs, double *u, double *v, int nsites, double *pr, size_t *ir,
          size_t *jc) {
  int k, iu, iv, col, Iu, Iv;
  int ret_code = 0;
  double nu[nsites];
  double nv[nsites];
  int interu[nsites];
  int interv[nsites];
  double basis_u[4];
  double basis_v[4];
  int npts = bbs->nptsu * bbs->nptsv;
  size_t nb_elem_col[npts];
  size_t cur_ind[npts];

  // Compute the normalized evaluation values and their interval numbers
  normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
  normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);

  // FILLING jc
  // First, compute the number of elements in each column
  // Check at the same time if the sites are inside the definition domain.
  for (k = 0; k < nsites; ++k) {
    Iu = interu[k];
    Iv = interv[k];

    if (Iu < 0 || Iu > bbs->nptsu - 4 || Iv < 0 || Iv > bbs->nptsv - 4) {
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

  for (k = 0; k < nsites; ++k) {
    eval_basis(nu[k], basis_u);
    eval_basis(nv[k], basis_v);
    for (iu = 0; iu < 4; ++iu)
      for (iv = 0; iv < 4; ++iv) {
        col = (iu + interu[k]) * bbs->nptsv + iv + interv[k];
        pr[cur_ind[col]] = basis_u[iu] * basis_v[iv];
        ir[cur_ind[col]] = k;
        cur_ind[col] += 1;
      }
  }

error:

  return ret_code;
}

/* Build a sparse colocation matrix that accounts for derivatives.
 * This function does not handle sites that are outside of the domain (an error
 * is issued in such cases).
 * ARGUMENTS:
 * - bbs: b-spline structure.
 * - u,v: sites used to build the colocation matrix.
 * - nsites: number of sites (number of elements in u and v).
 * - du,dv: order of derivation along the u- and the v- axis resp.
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
int coloc_deriv(bbs_t *bbs, double *u, double *v, int nsites, int du, int dv,
                double *pr, size_t *ir, size_t *jc) {
  int k, iu, iv, col, Iu, Iv;
  int ret_code = 0;
  double nu[nsites];
  double nv[nsites];
  int interu[nsites];
  int interv[nsites];
  basis_func_t b_func_u = get_basis_ptr(du);
  basis_func_t b_func_v = get_basis_ptr(dv);
  double fact = get_deriv_fact(bbs, du, dv);
  double basis_u[4];
  double basis_v[4];
  int npts = bbs->nptsu * bbs->nptsv;
  size_t nb_elem_col[npts];
  size_t cur_ind[npts];
  // Compute the normalized evaluation values and their interval numbers
  normalize_with_inter(bbs->umin, bbs->umax, bbs->nptsu, u, nsites, nu, interu);
  normalize_with_inter(bbs->vmin, bbs->vmax, bbs->nptsv, v, nsites, nv, interv);

  // FILLING jc
  // First, compute the number of elements in each column
  // Check at the same time if the sites are inside the definition domain.
  for (k = 0; k < nsites; ++k) {
    Iu = interu[k];
    Iv = interv[k];

    if (Iu < 0 || Iu > bbs->nptsu - 4 || Iv < 0 || Iv > bbs->nptsv - 4) {
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

  for (k = 0; k < nsites; ++k) {
    // eval_basis(nu[k], basis_u);
    // eval_basis(nv[k], basis_v);
    b_func_u(nu[k], basis_u);
    b_func_v(nv[k], basis_v);
    for (iu = 0; iu < 4; ++iu)
      for (iv = 0; iv < 4; ++iv) {
        col = (iu + interu[k]) * bbs->nptsv + iv + interv[k];
        pr[cur_ind[col]] = fact * basis_u[iu] * basis_v[iv];
        ir[cur_ind[col]] = k;
        cur_ind[col] += 1;
      }
  }

error:
  return ret_code;
}

/* Precomputed coefficients used in the bending matrix.
 * These arrays are used by the "bending_ur" function and they are
 * not intended to be used directly by the user. */
double __bxx_coeff[] = {
    1.0 / 756.0,     0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -1.0 / 504.0,    1.0 / 252.0,     0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             -1.0 / 504.0,    1.0 / 252.0,     0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    1.0 / 1512.0,    0.0,             -1.0 / 504.0,    1.0 / 756.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    43.0 / 5040.0,   -43.0 / 3360.0,  0.0,             43.0 / 10080.0,
    11.0 / 140.0,    0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -43.0 / 3360.0,  43.0 / 1680.0,   -43.0 / 3360.0,  0.0,
    -33.0 / 280.0,   33.0 / 140.0,    0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             -43.0 / 3360.0,  43.0 / 1680.0,   -43.0 / 3360.0,
    0.0,             -33.0 / 280.0,   33.0 / 140.0,    0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    43.0 / 10080.0,  0.0,             -43.0 / 3360.0,  43.0 / 5040.0,
    11.0 / 280.0,    0.0,             -33.0 / 280.0,   11.0 / 140.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    1.0 / 252.0,     -1.0 / 168.0,    0.0,             1.0 / 504.0,
    311.0 / 5040.0,  -311.0 / 3360.0, 0.0,             311.0 / 10080.0,
    11.0 / 140.0,    0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -1.0 / 168.0,    1.0 / 84.0,      -1.0 / 168.0,    0.0,
    -311.0 / 3360.0, 311.0 / 1680.0,  -311.0 / 3360.0, 0.0,
    -33.0 / 280.0,   33.0 / 140.0,    0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             -1.0 / 168.0,    1.0 / 84.0,      -1.0 / 168.0,
    0.0,             -311.0 / 3360.0, 311.0 / 1680.0,  -311.0 / 3360.0,
    0.0,             -33.0 / 280.0,   33.0 / 140.0,    0.0,
    0.0,             0.0,             0.0,             0.0,
    1.0 / 504.0,     0.0,             -1.0 / 168.0,    1.0 / 252.0,
    311.0 / 10080.0, 0.0,             -311.0 / 3360.0, 311.0 / 5040.0,
    11.0 / 280.0,    0.0,             -33.0 / 280.0,   11.0 / 140.0,
    0.0,             0.0,             0.0,             0.0,
    1.0 / 15120.0,   -1.0 / 10080.0,  0.0,             1.0 / 30240.0,
    1.0 / 252.0,     -1.0 / 168.0,    0.0,             1.0 / 504.0,
    43.0 / 5040.0,   -43.0 / 3360.0,  0.0,             43.0 / 10080.0,
    1.0 / 756.0,     0.0,             0.0,             0.0,
    -1.0 / 10080.0,  1.0 / 5040.0,    -1.0 / 10080.0,  0.0,
    -1.0 / 168.0,    1.0 / 84.0,      -1.0 / 168.0,    0.0,
    -43.0 / 3360.0,  43.0 / 1680.0,   -43.0 / 3360.0,  0.0,
    -1.0 / 504.0,    1.0 / 252.0,     0.0,             0.0,
    0.0,             -1.0 / 10080.0,  1.0 / 5040.0,    -1.0 / 10080.0,
    0.0,             -1.0 / 168.0,    1.0 / 84.0,      -1.0 / 168.0,
    0.0,             -43.0 / 3360.0,  43.0 / 1680.0,   -43.0 / 3360.0,
    0.0,             -1.0 / 504.0,    1.0 / 252.0,     0.0,
    1.0 / 30240.0,   0.0,             -1.0 / 10080.0,  1.0 / 15120.0,
    1.0 / 504.0,     0.0,             -1.0 / 168.0,    1.0 / 252.0,
    43.0 / 10080.0,  0.0,             -43.0 / 3360.0,  43.0 / 5040.0,
    1.0 / 1512.0,    0.0,             -1.0 / 504.0,    1.0 / 756.0};
double __bxy_coeff[] = {
    1.0 / 200.0,     0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    7.0 / 1200.0,    17.0 / 600.0,    0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -1.0 / 100.0,    -29.0 / 1200.0,  17.0 / 600.0,    0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -1.0 / 1200.0,   -1.0 / 100.0,    7.0 / 1200.0,    1.0 / 200.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    7.0 / 1200.0,    49.0 / 7200.0,   -7.0 / 600.0,    -7.0 / 7200.0,
    17.0 / 600.0,    0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    49.0 / 7200.0,   119.0 / 3600.0,  -203.0 / 7200.0, -7.0 / 600.0,
    119.0 / 3600.0,  289.0 / 1800.0,  0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -7.0 / 600.0,    -203.0 / 7200.0, 119.0 / 3600.0,  49.0 / 7200.0,
    -17.0 / 300.0,   -493.0 / 3600.0, 289.0 / 1800.0,  0.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -7.0 / 7200.0,   -7.0 / 600.0,    49.0 / 7200.0,   7.0 / 1200.0,
    -17.0 / 3600.0,  -17.0 / 300.0,   119.0 / 3600.0,  17.0 / 600.0,
    0.0,             0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -1.0 / 100.0,    -7.0 / 600.0,    1.0 / 50.0,      1.0 / 600.0,
    -29.0 / 1200.0,  -203.0 / 7200.0, 29.0 / 600.0,    29.0 / 7200.0,
    17.0 / 600.0,    0.0,             0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    -7.0 / 600.0,    -17.0 / 300.0,   29.0 / 600.0,    1.0 / 50.0,
    -203.0 / 7200.0, -493.0 / 3600.0, 841.0 / 7200.0,  29.0 / 600.0,
    119.0 / 3600.0,  289.0 / 1800.0,  0.0,             0.0,
    0.0,             0.0,             0.0,             0.0,
    1.0 / 50.0,      29.0 / 600.0,    -17.0 / 300.0,   -7.0 / 600.0,
    29.0 / 600.0,    841.0 / 7200.0,  -493.0 / 3600.0, -203.0 / 7200.0,
    -17.0 / 300.0,   -493.0 / 3600.0, 289.0 / 1800.0,  0.0,
    0.0,             0.0,             0.0,             0.0,
    1.0 / 600.0,     1.0 / 50.0,      -7.0 / 600.0,    -1.0 / 100.0,
    29.0 / 7200.0,   29.0 / 600.0,    -203.0 / 7200.0, -29.0 / 1200.0,
    -17.0 / 3600.0,  -17.0 / 300.0,   119.0 / 3600.0,  17.0 / 600.0,
    0.0,             0.0,             0.0,             0.0,
    -1.0 / 1200.0,   -7.0 / 7200.0,   1.0 / 600.0,     1.0 / 7200.0,
    -1.0 / 100.0,    -7.0 / 600.0,    1.0 / 50.0,      1.0 / 600.0,
    7.0 / 1200.0,    49.0 / 7200.0,   -7.0 / 600.0,    -7.0 / 7200.0,
    1.0 / 200.0,     0.0,             0.0,             0.0,
    -7.0 / 7200.0,   -17.0 / 3600.0,  29.0 / 7200.0,   1.0 / 600.0,
    -7.0 / 600.0,    -17.0 / 300.0,   29.0 / 600.0,    1.0 / 50.0,
    49.0 / 7200.0,   119.0 / 3600.0,  -203.0 / 7200.0, -7.0 / 600.0,
    7.0 / 1200.0,    17.0 / 600.0,    0.0,             0.0,
    1.0 / 600.0,     29.0 / 7200.0,   -17.0 / 3600.0,  -7.0 / 7200.0,
    1.0 / 50.0,      29.0 / 600.0,    -17.0 / 300.0,   -7.0 / 600.0,
    -7.0 / 600.0,    -203.0 / 7200.0, 119.0 / 3600.0,  49.0 / 7200.0,
    -1.0 / 100.0,    -29.0 / 1200.0,  17.0 / 600.0,    0.0,
    1.0 / 7200.0,    1.0 / 600.0,     -7.0 / 7200.0,   -1.0 / 1200.0,
    1.0 / 600.0,     1.0 / 50.0,      -7.0 / 600.0,    -1.0 / 100.0,
    -7.0 / 7200.0,   -7.0 / 600.0,    49.0 / 7200.0,   7.0 / 1200.0,
    -1.0 / 1200.0,   -1.0 / 100.0,    7.0 / 1200.0,    1.0 / 200.0};
double __byy_coeff[] = {
    1.0 / 756.0,    0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    43.0 / 5040.0,  11.0 / 140.0,    0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    1.0 / 252.0,    311.0 / 5040.0,  11.0 / 140.0,    0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    1.0 / 15120.0,  1.0 / 252.0,     43.0 / 5040.0,   1.0 / 756.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 504.0,   -43.0 / 3360.0,  -1.0 / 168.0,    -1.0 / 10080.0,
    1.0 / 252.0,    0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -43.0 / 3360.0, -33.0 / 280.0,   -311.0 / 3360.0, -1.0 / 168.0,
    43.0 / 1680.0,  33.0 / 140.0,    0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 168.0,   -311.0 / 3360.0, -33.0 / 280.0,   -43.0 / 3360.0,
    1.0 / 84.0,     311.0 / 1680.0,  33.0 / 140.0,    0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 10080.0, -1.0 / 168.0,    -43.0 / 3360.0,  -1.0 / 504.0,
    1.0 / 5040.0,   1.0 / 84.0,      43.0 / 1680.0,   1.0 / 252.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 504.0,   -43.0 / 3360.0,  -1.0 / 168.0,    -1.0 / 10080.0,
    1.0 / 252.0,    0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -43.0 / 3360.0, -33.0 / 280.0,   -311.0 / 3360.0, -1.0 / 168.0,
    43.0 / 1680.0,  33.0 / 140.0,    0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 168.0,   -311.0 / 3360.0, -33.0 / 280.0,   -43.0 / 3360.0,
    1.0 / 84.0,     311.0 / 1680.0,  33.0 / 140.0,    0.0,
    0.0,            0.0,             0.0,             0.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 10080.0, -1.0 / 168.0,    -43.0 / 3360.0,  -1.0 / 504.0,
    1.0 / 5040.0,   1.0 / 84.0,      43.0 / 1680.0,   1.0 / 252.0,
    0.0,            0.0,             0.0,             0.0,
    1.0 / 1512.0,   43.0 / 10080.0,  1.0 / 504.0,     1.0 / 30240.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 504.0,   -43.0 / 3360.0,  -1.0 / 168.0,    -1.0 / 10080.0,
    1.0 / 756.0,    0.0,             0.0,             0.0,
    43.0 / 10080.0, 11.0 / 280.0,    311.0 / 10080.0, 1.0 / 504.0,
    0.0,            0.0,             0.0,             0.0,
    -43.0 / 3360.0, -33.0 / 280.0,   -311.0 / 3360.0, -1.0 / 168.0,
    43.0 / 5040.0,  11.0 / 140.0,    0.0,             0.0,
    1.0 / 504.0,    311.0 / 10080.0, 11.0 / 280.0,    43.0 / 10080.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 168.0,   -311.0 / 3360.0, -33.0 / 280.0,   -43.0 / 3360.0,
    1.0 / 252.0,    311.0 / 5040.0,  11.0 / 140.0,    0.0,
    1.0 / 30240.0,  1.0 / 504.0,     43.0 / 10080.0,  1.0 / 1512.0,
    0.0,            0.0,             0.0,             0.0,
    -1.0 / 10080.0, -1.0 / 168.0,    -43.0 / 3360.0,  -1.0 / 504.0,
    1.0 / 15120.0,  1.0 / 252.0,     43.0 / 5040.0,   1.0 / 756.0};

/* This function computes the upper right part of the "bending matrix".
 * - bbs: the bidimensional b-spline structure.
 * - lambdas: array storing an (nptsv-3)x(nptsu-3) matrix by column
 *   lambdas[iu*(nptsv-3)+iv] is the regularization parameter over the knot
 * interval #(iu,iv).
 * - pr, ir, jc: internal storage for the resulting sparse bending matrix.
 */
void bending_ur(bbs_t *bbs, double err, double *pr, size_t *ir, size_t *jc) {
  int ny = bbs->nptsu; // Order of the dimensions are switched due to a previous
                       // version of the code
  int nx = bbs->nptsv;
  int i, j, a, b, c, d, e1, f1, e2, f2, curnb, total, ind, sb, nb;
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
  for (j = 0; j < ny; ++j) {
    for (i = 0; i < nx; ++i) {
      curnb = 0;
      for (b = max(j - 3, 0); b <= j - 1; ++b)
        for (a = max(0, i - 3); a <= min(nx - 1, i + 3); ++a) {
          *pi = b * nx + a;
          *px = 0.0;
          ++pi;
          ++px;
          ++curnb;
        }

      for (a = max(0, i - 3); a <= i; ++a) {
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
    for (i = 0; i <= j; ++i) {
      ind = 16 * j + i;
      coeff_b[ind] = sy * __bxx_coeff[ind] / pow(sx, 3) +
                     __bxy_coeff[ind] / (sx * sy) +
                     sx * __byy_coeff[ind] / pow(sy, 3);
    }

  // Put the right coefficients at the right locations (one knot domain at a
  // time
  // and with the scaling given by lambda)
  pp = jc;
  px = pr;
  for (b = 0; b < ny - 3; ++b) {
    for (a = 0; a < nx - 3; ++a) {
      lbd = err;
      for (c = 0; c < 16; ++c) {
        for (d = c; d < 16; ++d) {
          e1 = c / 4;
          f1 = c % 4;
          e2 = d / 4;
          f2 = d % 4;
          i = (b + e1) * nx + a + f1;
          j = (b + e2) * nx + a + f2;

          // See (notebook 2, p. 61) for the tricky formulas
          nb = i / nx - max(j / nx - 3, 0);
          sb = min(min(4 + (j % nx), 3 + nx - (j % nx)), min(nx, 7));
          px[pp[j] + nb * sb + (i % nx) - max((j % nx) - 3, 0)] +=
              lbd * coeff_b[16 * d + c];
        }
      }
    }
  }
}
} // End namespace
