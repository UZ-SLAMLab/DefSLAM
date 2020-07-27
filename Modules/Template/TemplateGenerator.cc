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

#include <TemplateGenerator.h>
#include <mutex>

#include "LaplacianMesh.h"

namespace defSLAM
{
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapPoint;

  class LaplacianMesh;
  // Load template with a keyframe that contains a surface. Used for shape-from-template.
  Template *
  TemplateGenerator::ChargeLaplacianMesh(std::vector<std::vector<double>> &vertex,
                                         std::vector<std::vector<int>> &index,
                                         Map *map)
  {
    LaplacianMesh *lapmesh = new LaplacianMesh(vertex, index, map);

    return lapmesh;

    return (static_cast<LaplacianMesh *>(nullptr));
  }

  // Load template with a keyframe that contains a surface. Used during in DefMap.
  Template *
  TemplateGenerator::LaplacianMeshCreate(std::set<MapPoint *> &mspMapPoints,
                                         Map *map, KeyFrame *kF)
  {
    LaplacianMesh *lapmesh = new LaplacianMesh(mspMapPoints, map, kF);

    return lapmesh;

    return (static_cast<LaplacianMesh *>(nullptr));
  }

} // namespace defSLAM
