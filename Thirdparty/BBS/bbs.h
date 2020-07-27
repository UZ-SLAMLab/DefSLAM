/* bbs.h
 *
 * General routines for cubic bidimensional b-splines.
 * 
 * History
 *  2009/??/??: First version
 *  2010/12/15: Translation to CPP
 *              Adding OpenMP support
 *  2018/12/18: Adaptation to c++ for defSLAM (Jose Lamarca)
 *
 * (c)2009-2010, Florent Brunet.
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

#ifndef __BBSS_H__
#define __BBSS_H__

#include <stdlib.h>
#include <math.h>
#include <omp.h>

namespace BBS
{
    /* The B-Spline structure */
    typedef struct _bbs_t
    {
        double umin;
        double umax;
        int nptsu;
        double vmin;
        double vmax;
        int nptsv;
        int valdim;
    } bbs_t;

    void normalize(double xmin, double xmax, int npts, double *x, int nb_x, double *nx);
    void normalize_with_inter(double xmin, double xmax, int npts, double *x, int nb_x, double *nx, int *inter);

    void eval_basis(double nx, double *val_basis);
    void eval_basis_d(double nx, double *val_basis);
    void eval_basis_dd(double nx, double *val_basis);

    void eval(bbs_t *bbs, double *ctrlpts, double *u, double *v, int nval, double *val, int du, int dv);

    int coloc(bbs_t *bbs, double *u, double *v, int nsites, double *pr, size_t *ir, size_t *jc);
    //int coloc(bbs_t *bbs, double *u, double *v, int nsites, std::vector<double> *pr, std::vector<size_t> *ir, std::vector<size_t> *jc);
    int coloc_deriv(bbs_t *bbs, double *u, double *v, int nsites, int du, int dv, double *pr, size_t *ir, size_t *jc);

    void bending_ur(bbs_t *bbs, double err, double *pr, size_t *ir, size_t *jc);
} // namespace BBS
#endif
