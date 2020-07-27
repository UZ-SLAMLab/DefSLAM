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

#include <SmootherMLS.h>

#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>
#include <pcl/common/distances.h>

#include <pcl/filters/radius_outlier_removal.h>

// Constructor.
SmootherMLS::SmootherMLS(int polynomialOrder_, double searchRadius_)
    : polynomialOrder_(polynomialOrder_), searchRadius_(searchRadius_){

                                          };

using namespace pcl;

// Constructor.
std::vector<int> SmootherMLS::smoothPointCloud(const SmootherMLS::pointcloud &pcVector)
{
    pcl::PointCloud<PointXYZ>::Ptr cloud(new pcl::PointCloud<PointXYZ>());

    for (uint i(0); i < pcVector.size(); i++)
    {
        PointXYZ pclpoint;
        pclpoint.x = pcVector[i][0];
        pclpoint.y = pcVector[i][1];
        pclpoint.z = pcVector[i][2];
        cloud->push_back(pclpoint);
    }

    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointNormal> mls_points;
    // Set parameters
    pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
    mls.setInputCloud(cloud);
    mls.setComputeNormals(true); // To have the corresponding Indices
    mls.setPolynomialOrder(polynomialOrder_);
    mls.setSearchMethod(tree);
    mls.setSearchRadius(searchRadius_);
    mls.process(mls_points);
    auto indices = (*(mls.getCorrespondingIndices())).indices;
    SmootherMLS::pointcloud pcVectorO;
    pcVectorO.reserve(mls_points.size());

    for (uint i(0); i < pcVector.size(); i++)
    {
        std::vector<float> pt;
        pt.reserve(3);
        pt.push_back(mls_points[i].x);
        pt.push_back(mls_points[i].y);
        pt.push_back(mls_points[i].z);
        pcVectorO.push_back(std::move(pt));
    }

    return indices;
}

std::vector<int> SmootherMLS::outlierRemovalRadius(const SmootherMLS::pointcloud &pcVector)
{
    pcl::PointCloud<PointXYZ>::Ptr cloud(new pcl::PointCloud<PointXYZ>());

    for (uint i(0); i < pcVector.size(); i++)
    {
        PointXYZ pclpoint;
        pclpoint.x = pcVector[i][0];
        pclpoint.y = pcVector[i][1];
        pclpoint.z = pcVector[i][2];
        cloud->push_back(pclpoint);
    }

    // Set parameters
    pcl::RadiusOutlierRemoval<PointXYZ> outrem;
    outrem.setInputCloud(cloud);
    outrem.setRadiusSearch(searchRadius_);
    outrem.setMinNeighborsInRadius(20);
    std::vector<int> indices;
    outrem.filter(indices);
    return indices;
}