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

#ifndef TRIANGULARMESH_H
#define TRIANGULARMESH_H
#include "opencv2/opencv.hpp"
#include <Facet.h>
#include <KeyFrame.h>
#include <MapPoint.h>
#include <Optimizer.h>

#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>

#include "Template.h"
#include <Converter.h>

namespace defSLAM
{

  class Template;
  class MapPoint;
  class KeyFrame;
  class Frame;

  class TriangularMesh : public Template
  {
  public:
    // Constructor by default.
    TriangularMesh() = delete;

    // Constructor from mesh given nodes and facets.
    TriangularMesh(std::vector<std::vector<double>> &,
                   std::vector<std::vector<int>> &, Map *map);

    // Constructor for a mesh given a surface estimated in the kf.
    TriangularMesh(std::set<MapPoint *> &mspMapPoints, Map *map, KeyFrame *kf);

    // Destructor
    virtual ~TriangularMesh() = default;

  public:
    // Get the texture of a facet for drawing it in the viewer.
    void getFacetTexture(KeyFrame *KF);

  private:
    // Create a regular triangular mesh from the keyframe surface estimated.
    std::vector<std::vector<int>> regularTriangulation(int nodesVer, int nodesHor);

    // Create the nodes and the facets to define the triagular mesh.
    void setNodes(std::vector<std::vector<float>> &,
                  std::vector<std::vector<int>> &);

    // Embed the points inside the facets by calculating its barycentric coordinates.
    void calculateFeaturesCoordinates();
    //Estimate barycentric coordinates. It returns true if the point is inside the face.
    bool pointInTriangle(const Eigen::Vector3f &query_point,
                         const Eigen::Vector3f &triangle_vertex_0,
                         const Eigen::Vector3f &triangle_vertex_1,
                         const Eigen::Vector3f &triangle_vertex_2,
                         Eigen::Vector3f &barycentric);

    void discardFaces();
  };
} // namespace defSLAM
#endif
