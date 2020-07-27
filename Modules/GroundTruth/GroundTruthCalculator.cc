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

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "GroundTruthCalculator.h"
#include <Eigen/Dense>
#include <fstream>
#include "opencv2/calib3d.hpp"
#include "CC_MAC.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;

namespace defSLAM
{
  namespace GroundTruthTools
  {
    /**********************
 * The funtion ScaleMinMedian return the scale estimated
 * for two point clouds with a min median filter
 * 
 * Arg:
 *  PosMono std::vector<std::vector<floats>> points to estimated by DefSLAM
 *  PosStereo std::vector<std::vector<floats>> points estimated by stereo
 * 
 * Return:
 *  Best scale: best estimation of the scale
 * 
 * Author: Javier Morlana
 * 
 * *******************/
    float scaleMinMedian(
        std::vector<std::vector<float>> &PosMono,
        std::vector<std::vector<float>> &PosStereo)
    {

      float min_med = 10000.0;
      int n_points = PosMono.size();
      int final_points = 0;
      double best_scale(0.0);

      for (int i = 0; i < n_points; i++)
      {
        double r_i = ((double)rand() / (RAND_MAX));

        if (r_i > 0.25)
          continue;
        double scale =
            PosStereo[i][2] / PosMono[i][2]; // stereo_z/mono_z

        vector<float> squared_res;
        squared_res.resize(n_points);
        {
          for (int j = 0; j < n_points; j++)
          {
            squared_res[j] = -1;
            if (i == j)
              continue;
            double r_j = ((double)rand() / (RAND_MAX));

            if (r_j > 0.25)
              continue;
            float r2(0.0);
            for (uint k(0); k < 3; k++)
            {
              auto res = (scale * PosMono[j][k] - PosStereo[j][k]);
              r2 = r2 + res * res;
            }
            squared_res[j] = std::sqrt(r2);
          }
        }

        // Get the median of the residual
        std::sort(squared_res.begin(), squared_res.end());

        int NumberNonZero(0);

        while (squared_res[NumberNonZero++] < 0)
          continue;

        std::vector<float> NonZeroVec(squared_res.begin() + NumberNonZero,
                                      squared_res.end());

        int median_index = round(NonZeroVec.size() / 2);
        final_points++;

        if (NonZeroVec.size() == 0)
          return 0.0;

        if (NonZeroVec[median_index] < min_med)
        {
          // Save the minimum median
          min_med = NonZeroVec[median_index];
          best_scale = scale;
        }

        squared_res.clear();
      }

      // Desviation
      float desv = 1.4826 * (1.0 - (5.0 / (final_points - 1.0))) * sqrt(min_med);
      std::cout << " stan dev " << desv << std::endl;

      // Comparison and rejecting outliers
      std::vector<int> NonOutlier;
      for (int i = 0; i < n_points; i++)
      {
        float residual(0.0);
        for (uint j(0); j < 3; j++)
        {
          auto rest = (best_scale * PosMono[i][j] - PosStereo[i][j]);
          residual = residual + rest * rest;
        }
        residual = std::sqrt(residual);

        if ((residual / desv) < 2.5)
        {
          NonOutlier.push_back(i);
        }
      }

      n_points = NonOutlier.size();

      // Final scale with inliers
      float Sum_num = 0.0;
      float Sum_den = 0.0;

      for (int i = 0; i < n_points; i++)
      {
        int indx = NonOutlier[i];
        Sum_num += (PosStereo[indx][2] * PosMono[indx][2]);
        Sum_den += (PosMono[indx][2] * PosMono[indx][2]);
      }

      best_scale = Sum_num / Sum_den;
      return best_scale;
    }

    /**********************
 * The funtion SaveResults saves a vector in a text file
 * 
 * Author: Jose Lamarca
 * 
 * Arg:
 *  vectorError std::vector<floats> save a vector 
 *  name std::string name for the file
 * 
 * Return:
 *  Best scale: best estimation of the scale
 * 
 * *******************/
    void saveResults(std::vector<float> &vectorError,
                     std::string &name)
    {
      Eigen::MatrixXd asd(vectorError.size(), 1);
      for (uint i(0); i < vectorError.size(); i++)
      {
        asd(i, 0) = vectorError[i];
      }
      std::ofstream myfile;
      myfile.open(name);
      myfile << asd;
      myfile.close();
    }
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
 *  name std::string name for the file
 * 
 * Return:
 *  std::vector<float> point in 3D in camera frame;
 * 
 * *******************/
    std::vector<float> estimateGT(const cv::KeyPoint &kp, const cv::Mat &imLeft, const cv::Mat &imRight, double mbf,
                                  double cx, double cy, double fx, double fy)
    {
      if (kp.pt.x < 0)
        return std::vector<float>(3, -1);
      if (kp.pt.y < 0)
        return std::vector<float>(3, -1);
      if (kp.pt.y > imRight.rows - 0)
        return std::vector<float>(3, -1);
      if (kp.pt.x > imRight.cols - 60)
        return std::vector<float>(3, -1);
      int tempx(TEMPX), tempy(TEMPY);
      int searchx(SEARCHX), searchy(tempy + MARGIN);

      if (((kp.pt.x - tempx / 2) < 20) or ((kp.pt.y - tempy / 2) < 0) or
          ((kp.pt.x + tempx / 2) > imRight.cols) or
          ((kp.pt.y + tempy / 2) > imRight.rows))
        return std::vector<float>(3, -1);

      cv::Rect cropRect(kp.pt.x - tempx / 2, kp.pt.y - tempy / 2, tempx, tempy);

      int finx = (kp.pt.x - searchx);
      int finy = (kp.pt.y - searchy);
      int finxright = kp.pt.x + mbf / 4;
      int finyright = finy + searchy * 2;
      if (finx < 0)
        finx = 0;
      if (finy < 0)
        finy = 0;
      if (finxright > imRight.cols)
        searchx = float(imRight.cols - 1 - finx);
      if (finyright > imRight.rows)
        searchy = float(imRight.rows - 1 - finy) / 2;

      cv::Rect SearchRegion(finx, finy, searchx, searchy * 2);
      cv::Mat tmp = imLeft(cropRect);
      double minValcrop;
      double maxValcrop;
      cv::Point minLoccrop;
      cv::Point maxLoccrop;
      cv::minMaxLoc(tmp, &minValcrop, &maxValcrop, &minLoccrop, &maxLoccrop, cv::Mat());
      if (maxValcrop > 250)
        return std::vector<float>(3, -1);
      cv::Mat SearchImage = imRight(SearchRegion);
      cv::Mat result;
      cv::matchTemplate(SearchImage, tmp, result, cv::TM_CCORR_NORMED);
      cv::Point matchLoc;

      double minVal;
      double maxVal;
      cv::Point minLoc;
      cv::Point maxLoc;
      cv::minMaxLoc(result, &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat());

      double CorrelationThreshold(THRESHOLD);

      if (maxVal < CorrelationThreshold)
        return std::vector<float>(3, -1);

      matchLoc = maxLoc;

      matchLoc.x = matchLoc.x + finx + tempx / 2;
      matchLoc.y = matchLoc.y + finy + tempy / 2;

      float disp(std::abs(matchLoc.x - kp.pt.x));
      std::vector<float> ps;
      ps.reserve(3);
      ps.push_back(mbf / disp * (((float)kp.pt.x - cx) / fx));
      ps.push_back(mbf / disp * (((float)kp.pt.y - cy) / fy));
      ps.push_back(mbf / disp);
      return ps;
    }
  } // namespace GroundTruthTools
} // namespace defSLAM
