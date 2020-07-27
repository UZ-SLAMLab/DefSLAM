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

#ifndef TEMPLATEGENERATOR_H
#define TEMPLATEGENERATOR_H

#include "KeyFrame.h"
#include "Map.h"
#include "MapPoint.h"
#include "Template.h"

namespace ORB_SLAM2
{
  class MapPoint;
  class KeyFrame;
} // namespace ORB_SLAM2
namespace defSLAM
{
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::MapPoint;

  class TemplateGenerator
  {
  public:
    //// Load a template given the nodes and facets.
    static Template *ChargeLaplacianMesh(std::vector<std::vector<double>> &,
                                         std::vector<std::vector<int>> &,
                                         Map *map);

    // Load template with a keyframe that contains a surface.
    static Template *LaplacianMeshCreate(std::set<MapPoint *> &mspMapPoints,
                                         Map *map, KeyFrame *kF);
  };

} // namespace defSLAM

#endif // OPTIMIZER_H
