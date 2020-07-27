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
#ifndef POLYSOLVER_H
#define POLYSOLVER_H
#include <Eigen/Core>
#include <ceres/ceres.h>
#include <iostream>

namespace defSLAM
{

  struct Eqs
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Matrix<double, 1, 10> eqs1;
    Eigen::Matrix<double, 1, 10> eqs2;
  };

  class PolySolver : public ceres::SizedCostFunction<2, 2>
  {
    /***********
    * This class is a block of ceres used to estimate the NRSfM.
    * The equations used in this class use the warp estimated from
    * the image one to the image 2 \eta_{12}. There is jupyter notebook
    * and a script in python to develop the polynomials with a simbolic 
    * library. 
   * *******/

  public:
    // Initialize
    PolySolver(const Eigen::Matrix<double, 1, 10> &eq1,
               const Eigen::Matrix<double, 1, 10> &eq2)
    {
      init_d(eq1, eq2);
    }
    ~PolySolver() { delete eqs; }

    // Evaluate function for Ceres
    bool Evaluate(double const *const *parameters, double *e,
                  double **Jac) const;

    /****
     * This function gets the differential parameters and 
     * returns the coefficients of the polynomials. There 
     * are two polinomials. Eq (13) and (14) of the paper.
     * *****/
    static void getCoefficients(double a, double b, double c, double d,
                                double t1, double t2, double e, double e1,
                                double x1, double y1, double x2, double y2,
                                int i, double *eq);

  private:
    // Helper to initialize Eigen parameters in a class. It gave
    // some problems.
    void init_d(const Eigen::Matrix<double, 1, 10> &eq1,
                const Eigen::Matrix<double, 1, 10> &eq2)
    {
      eqs = new Eqs;
      eqs->eqs1 = eq1;
      eqs->eqs2 = eq2;
    }
    Eqs *eqs;
  };
} // namespace defSLAM
#endif
