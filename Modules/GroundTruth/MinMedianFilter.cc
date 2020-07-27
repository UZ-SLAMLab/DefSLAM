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
#include <cmath>
#include <algorithm>
#include <stdio.h>

#include "MinMedianFilter.h"

using namespace std;

namespace defSLAM
{
  namespace PointCloudFilters
  {
    /**********************
 * The funtion minMedianFilter return the non-outliers in a ditribution
 * applying a min median filter. 
 * Arg:
 *  PointsToFilter std::vector<std::vector<floats>> points to filter
 * Return:
 *  NonOutlier: Index of the points in the distribution
 * 
 * *******************/
    vector<int> minMedianFilter(vector<vector<float>> PointsToFilter)
    {
      float min_med = 10000.0;
      int n_points = PointsToFilter.size();
      int final_points = 0;
      double best_depth(0.0);
      std::vector<int> NonOutlier;
      NonOutlier.reserve(n_points);
      for (int i = 0; i < n_points; i++)
      {
        vector<float> squared_res;

        double r_i = ((double)rand() / (RAND_MAX));
        if (r_i > 0.75)
          continue;
        float depth = PointsToFilter[i][2];

        if (depth == -1.0)
          continue;

        for (int j = 0; j < n_points; j++)
        {
          if (i == j)
            continue;
          double r_j = ((double)rand() / (RAND_MAX));
          if (r_j > 0.75)
            continue;
          float depth_2 = PointsToFilter[j][2];
          if (depth_2 == -1.0)
            continue;
          float r2 = (depth - depth_2);
          squared_res.push_back(r2 * r2);
        }
        // Get the median of the residual
        sort(squared_res.begin(), squared_res.end());
        int median_index = round(squared_res.size() / 2);

        final_points++;

        if (squared_res.size() == 0)
          return NonOutlier;

        if (squared_res[median_index] < min_med)
        {
          // Save the minimum median
          min_med = squared_res[median_index];
          best_depth = depth;
        }
        squared_res.clear();
      }

      std::cout << " final points and min median " << final_points << " " << min_med
                << std::endl;

      // Desviation
      float desv = 1.4826 * (1.0 - (5.0 / (final_points - 1.0))) * std::sqrt(min_med);
      // Comparison and rejecting outliers
      for (int i = 0; i < n_points; i++)
      {
        float depth = PointsToFilter[i][2];
        if (depth == -1.0)
          continue;

        float residual = best_depth - depth;

        residual = residual * residual;
        if ((residual / desv) < 2.0)
        {
          NonOutlier.push_back(i);
        }
      }
      n_points = NonOutlier.size();

      //  cout << "MIN MEDIAN FILTER: inliers " << n_points << endl;
      return NonOutlier;
    }
  } // namespace PointCloudFilters
} // namespace defSLAM
