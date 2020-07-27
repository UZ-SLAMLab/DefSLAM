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

#ifndef BBSCOLOC_H
#define BBSCOLOC_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <bbs.h>

namespace BBS
{
    void colocEigen(bbs_t *bbs, double *u, double *v, int nsites,
                    Eigen::MatrixXd &colocMatrix);

    int coloc_derivEigen(bbs_t *bbs, double *u, double *v, int nsites, int du,
                         int dv, Eigen::MatrixXd &colocMatrix);

    int coloc_derivEigen(bbs_t *bbs, double *u, double *v, int nsites, int du,
                         int dv, Eigen::MatrixXf &colocMatrix);

    void BendingEigen(bbs_t *bbs, double err,
                      Eigen::SparseMatrix<double> &benMatrix);

    void BendingEigen(bbs_t *bbs, double err, Eigen::MatrixXd &benMatrix);

    void EvalEigen(bbs_t *bbs, const double *ctrlpts, double *u, double *v,
                   int nval, Eigen::MatrixXd &val, int du, int dv);
} // namespace BBS

#endif
