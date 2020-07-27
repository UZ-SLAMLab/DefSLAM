
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

#ifndef LAPLACIANMESH_H
#define LAPLACIANMESH_H

#include "KeyFrame.h"
#include "MapPoint.h"
#include "TriangularMesh.h"
#include "set"

namespace Eigen
{
  typedef Matrix<double, 1, 1> Vector1d;
}

namespace ORB_SLAM2
{
  class TriangularMesh;
  class MapPoint;
  class KeyFrame;
} // namespace ORB_SLAM2
namespace defSLAM
{
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapPoint;

  class LaplacianMesh : public TriangularMesh
  {
  public:
    // Constructor from mesh.
    LaplacianMesh(std::vector<std::vector<double>> &vertex,
                  std::vector<std::vector<int>> &index, Map *map);

    // Constructor from surface estimated for kF.
    LaplacianMesh(std::set<MapPoint *> &mspMapPoints, Map *map, KeyFrame *kF);

    // Destructor
    virtual ~LaplacianMesh();

    /********
   * This function extract mean curvatures. It is thought for irregular meshes.
   * it use the method proposed in to extract Laplacian coordiantes that represent
   * the mean curvature.
   * ********/
    void ExtractMeanCurvatures();

    // Get the initial laplacian coordinate(vector of mean curvature) of the node
    const Eigen::Vector3d GetLaplacianCoord(Node *);

    // Get the initial mean curvature (norm of mean curvature)) of the node
    Eigen::Vector1d const GetMeanCurvatureInitial(Node *n);

    // Get current mean curvature. It consider that the weights remain constant.
    const Eigen::Vector1d GetMeanCurvature(Node *n);

    // Get current mean curvature as Laplacian vector
    Eigen::Vector3d const GetMeanCurvatureVector(Node *n);

    std::map<Node *, Eigen::Vector3d, std::less<Node *>,
             Eigen::aligned_allocator<std::pair<const Node *, Eigen::Vector3d>>>
        LaplacianCoords; // Laplacian coordinates.
  };
} // namespace defSLAM
#endif
