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
#ifndef MINMEDIANFILTER_H
#define MINMEDIANFILTER_H

#include <stdio.h>
#include <iostream>
#include <vector>

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
 *  NonOutlier: Index of the points that are in the distribution
 * 
 * *******************/
    std::vector<int> minMedianFilter(std::vector<std::vector<float>> PointsToFilter);
  } // namespace PointCloudFilters
} // namespace defSLAM

#endif // MINMEDIANFILTER_H
