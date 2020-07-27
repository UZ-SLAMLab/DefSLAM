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
#pragma once

#include <vector>
#include <iostream>

class SmootherMLS
{
    /************************************* 
     * Class to smooth using Moving Least Squares or outlierRemovalRadius.
     * 
     * args: 
     *    polynomialOrder_: polynomial fitting order
     *    searchRadius_: radius of point to fit the polinomial.
     * methods:
     *    smoothPointCloud:  function to smooth pointclouds
     *    input : std::vector<std::vector<float>> point cloud to smooth
     *    output : std::vector<std::vector<float>> smoothed point cloud 
     * 
     * Author: Jose Lamarca
     * Inspired in https://pcl-tutorials.readthedocs.io/en/latest/resampling.html?highlight=resampling#smoothing-and-normal-estimation-based-on-polynomial-reconstruction
     * 
     **********************************/
public:
    typedef std::vector<std::vector<float>> pointcloud;

    // Initialize the smoother.
    SmootherMLS(int polynomialOrder_ = 2, double searchRadius_ = 0.03);

    // Smooth a point cloud with MLS
    std::vector<int> smoothPointCloud(const pointcloud &);

    // Outlier Removal with radius criteria.
    std::vector<int> outlierRemovalRadius(const SmootherMLS::pointcloud &);

private:
    int polynomialOrder_; // polynomial order for the MLS
    double searchRadius_; // Seach radius.
};