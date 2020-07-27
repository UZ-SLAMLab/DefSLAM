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

#ifndef GROUNDTRUTHCALCULATOR_H
#define GROUNDTRUTHCALCULATOR_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include "opencv2/core/core.hpp"

namespace defSLAM
{
     namespace GroundTruthTools
     {
          /**********************
 * The funtion ScaleMinMedian return the scale estimated
 * for two point clouds with a min median filter
 * 
 * Author: Javier Morlana
 * 
 * Arg:
 *  PosMono std::vector<std::vector<floats>> points to estimated by DefSLAM
 *  PosStereo std::vector<std::vector<floats>> points estimated by stereo
 * 
 * Return:
 *  Best scale: best estimation of the scale
 * 
 * *******************/
          float scaleMinMedian(std::vector<std::vector<float>> &PosMono,
                               std::vector<std::vector<float>> &PosStereo);
          /**********************
 * The funtion SaveResults saves a vector in a text file
 * 
 * Author: Javier Morlana
 * 
 * Arg:
 *  vectorError std::vector<floats> save a vector 
 *  name std::string name for the file
 * 
 * Return:
 *  Best scale: best estimation of the scale
 * 
 * *******************/
          void saveResults(std::vector<float> &vectorError, std::string &name);
          /**********************
 * The funtion estimateGT estimates the 3D point from
 * a stereo pair and a keypoint
 * 
 * Author: Jose Lamarca
 * 
 * Arg:
 *  std::vector<float> kp2D
 *  im1 image of keypoint
 *  im2 image stereo to query
 *  mbf: Baseline of the stereo pair
 *  fx,fy,cx,cy: Focals and optical centers;
 * 
 * Return:
 *  std::vector<float> point in 3D in camera frame;
 * 
 * *******************/
          std::vector<float> estimateGT(const cv::KeyPoint &kp2D, const cv::Mat &im1, const cv::Mat &im2, double mbf,
                                        double cx, double cy, double fx, double fy);
     }; // namespace GroundTruthTools

} // namespace defSLAM

#endif // GROUNDTRUTHCALCULATOR_H
