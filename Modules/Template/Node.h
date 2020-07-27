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

#ifndef NODE_H
#define NODE_H

#include <Template.h>
#include <Edge.h>
#include <Facet.h>
#include <MapPoint.h>
#include <set>
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
  class Template;
  class Edge;
  class Facet;

  using ORB_SLAM2::Frame;
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapPoint;

  class Node : public MapPoint
  {
  public:
    // Constructor by default
    Node() = default;

    // Constructor given x,y,z, the indx and the template that contain
    // the nodes.
    Node(double x, double y, double z, uint indx, Template *);

    // Destructor
    ~Node();

    // Get index for the optimizer.
    uint getIndex();

    // Set index for the optimizer.
    void setIndex(uint i);

    // Add edge that contain this node.
    void addEdge(Edge *);

    // Remove edge that contain this node.
    void eraseEdge(Edge *edge);

    // Add facet that contain this node.
    void addFacet(Facet *);

    // Remove facet that contain this node.
    void eraseFacet(Facet *);

    // get edges that contain this node.
    std::set<Edge *> getEdges();

    // Get Facets connected with this node.
    std::set<Facet *> GetFacets();

    // Return distance from this node to another node
    double distanceto(Node *node);

    // Check that there is at least one facet that contain
    // this node.
    bool noFacets();

    // Set this node to remove.
    void setBadFlag();

    // Get nodes connected by an edge.
    std::set<Node *> GetNeighbours();

    // Set this node is in the Viewed zone.
    void setViewed();

    // Set that this node is in the Local zone.
    void setLocal();

    // Update bools viewed and local.
    void update();

    // Is this node in the viewed zone.
    bool isViewed();

    // Is this node in the local zone.
    bool isLocal();

    // Set role to NonObserved.
    void resetRole();

    // Reset node to shape at rest
    void reset();

    // Set 3D position of the node.
    void setXYZ(double xx, double yy, double zz);

    // Get 3D position of the node.
    void getXYZ(double &xx, double &yy, double &zz);

    // Get its initial 3D position.
    void getInitialPose(double &xx, double &yy, double &zz);

    // Set if the node is part of the boundary of the mesh.
    void setBoundary();

    // Is it boundary.
    bool isBoundary();

  public:
    std::map<Node *, double> weights;                       // Laplacian weights.
    std::map<Node *, std::pair<Node *, Node *>> NodesjJ_1J; // Ring of those nodes.
    double x, y, z;                                         // 3D pose.
    double xO, yO, zO;                                      // 3D initial pose.
    float proju, projv;                                     // Projection in the keyframe image.
    int posVecDrawer;                                       // Position in the mesh drawer.

  private:
    uint indx;                // Index in the optimization.
    Template *mTemplate;      // Template.
    std::set<Edge *> mEdge;   // Edges containing this node.
    std::set<Facet *> mFacet; // Facets containing this node.

    bool viewed, local; // Is this node viewed? and local?
    enum eROLE
    {
      VIEWED = 0, // Viewed zone implied that is in a facet with observations.
      LOCAL = 1,  // Local zone implied that is in a facet neighbour of facets with observations.
      NONOBS = 2
    };

    eROLE role;
    bool Boundary; // Is part of the boundary of the mesh.
    std::mutex mMutexPos;
    std::mutex mMutexCond;
  };

} // namespace defSLAM

#endif // NODE_H
