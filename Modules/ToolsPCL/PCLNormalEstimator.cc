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
#include "PCLNormalEstimator.h"

#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>

// Constructor from pcl to estimate the normals. It directly estimate the normals.
// The radious search the neighbours in the point cloud to estimate the normals.
// Points introduced as [X0,Y0,Z0,X1,Y1,Z1,...,Xi,Yi,Zi] in the same vector
PCLNormalEstimator::PCLNormalEstimator(std::vector<float> &Points, double radious)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);

    for (size_t i = 0; i < Points.size(); i = i + 3)
    {
        cloud->push_back(pcl::PointXYZ(Points[i], Points[i + 1], Points[i + 2]));
    }

    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    ne.setInputCloud(cloud);

    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
    ne.setSearchMethod(tree);
    ne.setRadiusSearch(radious);
    ne.setViewPoint(0, 0, 0);

    ne.compute(*cloud_normals);

    for (uint i(0); i < cloud_normals->size(); i++)
    {
        normals.push_back(cloud_normals->at(i).normal_x);
        normals.push_back(cloud_normals->at(i).normal_y);
        normals.push_back(cloud_normals->at(i).normal_z);
    }
}

// Get the result.
// Normals given as [X0,Y0,Z0,X1,Y1,Z1,...,Xi,Yi,Zi] in the same vector
std::vector<float> &&PCLNormalEstimator::getNormals()
{
    return std::move(normals);
}
