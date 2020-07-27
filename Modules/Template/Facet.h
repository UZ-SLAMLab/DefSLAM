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

#ifndef FACET_H
#define FACET_H

#include <Template.h>
#include <KeyFrame.h>
#include <Edge.h>
#include <Node.h>
#include <MapPoint.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>

#include <mutex>

namespace ORB_SLAM2
{
  class KeyFrame;
  class Map;
  class Frame;
  class MapPoint;
} // namespace ORB_SLAM2
namespace defSLAM
{
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapPoint;

  class Template;
  class Node;
  class Edge;

  class Facet
  {
  public:
    // Constructor. It create a facet from three ptr to nodes for the
    // template temp.
    Facet(Node *, Node *, Node *, Template *);

    // Destructor
    ~Facet();

    // get nodes in a facet as a set.
    std::set<Node *> getNodes();

    // get nodes in a facet as an array.
    std::array<Node *, 3> getNodesArray();

    // Check there are observations in the facet.
    bool noObservations();

    // Add map point observed in the facet.
    void addMapPoint(MapPoint *);

    // Erase map point observed in the facet.
    void eraseMapPoint(MapPoint *);

    // Get edges of the facet.
    std::set<defSLAM::Edge *> getEdges();

    // Get map points in the facet.
    std::set<MapPoint *> getMapPoints();

    // Remove facet
    void setBadFlag();

    // Get texture for the image in the keyframe given
    void getTextureCoordinates(KeyFrame *KF);

  public:
    std::pair<KeyFrame *, std::map<Node *, cv::Point>> texture;     // Image coordinates for the texture in the viewer.
    double R, G, B;                                                 // In case of no texture, it is drawed in gray.
    static std::map<KeyFrame *, std::set<Facet *>> TexturesDataset; // texture
                                                                    // for the facets.

  private:
    Template *template_;          // template this face belongs to.
    std::set<Node *> Nodes;       // Nodes in  this face as set.
    std::array<Node *, 3> mNodes; // Nodes in this face as array
    std::set<Edge *> edges_;      // Edges in this face as set.
    std::set<MapPoint *> mPoints; // Map points contained in this facet
    std::mutex mMutexKeyframe;
  };

} // namespace defSLAM

#endif // FACET_H
