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

#include "Node.h"

#include <mutex>

namespace defSLAM
{
  // Constructor.
  Node::Node(double x, double y, double z, uint indx, Template *temp)
      : x(x), y(y), z(z), xO(x), yO(y), zO(z), proju(0), projv(0), posVecDrawer(-1), indx(indx), mTemplate(temp),
        role(NONOBS), Boundary(false)
  {
    viewed = false;
    local = false;
  }

  // Destructor.
  Node::~Node()
  {
    mTemplate = nullptr;
    weights.clear();
    NodesjJ_1J.clear();
  }

  // Get index in the optimization. This index
  // makes things easier to work with the optimizer.
  uint Node::getIndex() { return indx; }

  // Set index in the optimization.
  void Node::setIndex(uint i) { indx = i; }

  // Add edge that contain this node.
  void Node::addEdge(Edge *edge) { mEdge.insert(edge); }

  // get edges that contain this node.
  void Node::eraseEdge(Edge *edge) { mEdge.erase(edge); }

  // Add facet that contain this node.
  void Node::addFacet(Facet *facet) { mFacet.insert(facet); }

  // Remove facet that contain this node.
  void Node::eraseFacet(Facet *facet) { mFacet.erase(facet); }

  // get edges that contain this node.
  std::set<Edge *> Node::getEdges() { return mEdge; }

  // Get Facets connected with this node.
  std::set<Facet *> Node::GetFacets() { return mFacet; }

  // Return distance from this node to another node
  double Node::distanceto(Node *node)
  {
    double d = pow(this->x - node->x, 2) + pow(this->y - node->y, 2) +
               pow(this->z - node->z, 2);
    return sqrt(d);
  }

  // Check that there is at least one facet that contain
  // this node.
  bool Node::noFacets() { return (mFacet.size() == 0); }

  // Update the role
  void Node::update()
  {
    std::unique_lock<std::mutex> a(mMutexCond);
    if (VIEWED == role)
    {
      local = false;
      viewed = true;
    }
    else if (LOCAL == role)
    {
      viewed = false;
      local = true;
    }
    else
    {
      viewed = false;
      local = false;
    }
  }

  // Set this node to remove.
  void Node::setBadFlag()
  {
    mTemplate->eraseNode(this);
    for (std::set<Edge *>::iterator it = mEdge.begin(); it != mEdge.end(); it++)
    {
      (*it)->setBadFlag();
    }
  }

  // Get nodes connected by an edge.
  std::set<Node *> Node::GetNeighbours()
  {
    std::set<Node *> Neighbours;
    for (std::set<Edge *>::iterator Edges = mEdge.begin(); Edges != mEdge.end();
         Edges++)
    {
      std::set<Node *> EdgeNodes = (*Edges)->getNodes();
      for (std::set<Node *>::iterator Nodes = EdgeNodes.begin();
           Nodes != EdgeNodes.end(); Nodes++)
      {
        if (*Nodes != this)
          Neighbours.insert(*Nodes);
      }
    }
    return Neighbours;
  }

  // Set this node is in the Viewed zone.
  void Node::setViewed()
  {
    std::unique_lock<std::mutex> a(mMutexCond);
    role = VIEWED;
  }

  void Node::setLocal()
  {
    std::unique_lock<std::mutex> a(mMutexCond);
    if (role != VIEWED)
      role = LOCAL;
  }

  void Node::resetRole()
  {
    std::unique_lock<std::mutex> a(mMutexCond);
    role = NONOBS;
  }

  // Reset node to shape at rest
  void Node::reset()
  {
    std::unique_lock<mutex> as(mMutexPos);
    x = xO;
    y = yO;
    z = zO;
  }

  // Is this node in the viewed zone.
  bool Node::isViewed()
  {
    std::unique_lock<std::mutex> a(mMutexCond);
    return viewed;
  }

  // Is this node in the local zone.
  bool Node::isLocal()
  {
    std::unique_lock<std::mutex> a(mMutexCond);
    return local;
  }

  // Set 3D position of the node.
  void Node::setXYZ(double xx, double yy, double zz)
  {
    std::unique_lock<std::mutex> a(mMutexPos);
    this->x = xx;
    this->y = yy;
    this->z = zz;
  }

  // Get 3D position of the node.
  void Node::getXYZ(double &xx, double &yy, double &zz)
  {
    std::unique_lock<std::mutex> a(mMutexPos);
    xx = this->x;
    yy = this->y;
    zz = this->z;
  }

  // Get its initial 3D position.
  void Node::getInitialPose(double &xx, double &yy, double &zz)
  {
    xx = xO;
    yy = yO;
    zz = zO;
  }

  // Set if the node is part of the boundary of the mesh.
  void Node::setBoundary() { Boundary = true; }

  // Is the node in the boundary.
  bool Node::isBoundary() { return Boundary; }

} // namespace defSLAM
