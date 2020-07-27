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

#include <Edge.h>
#include <mutex>

namespace defSLAM
{
  // Constructor. Initialize an edge from two nodes, a facet
  // and the template
  Edge::Edge(Node *v1, Node *v2, Facet *facet, Template *temp)
  {
    mNodes[0] = v1;
    mNodes[1] = v2;
    mspNode.insert(v1);
    mspNode.insert(v2);

    bool RepeatedEdge(false);
    const std::set<Edge *> &MapEdges = temp->getEdges();
    // Check if the edge exist in the new template
    for (std::set<Edge *>::const_iterator it = MapEdges.begin();
         it != MapEdges.end(); it++)
    {
      if ((*it)->operator==(*this))
      {
        RepeatedEdge = true;
        break;
      }
    }

    //  If the edge does not exist, it is created.
    if (!RepeatedEdge)
    {
      InitialDist = v1->distanceto(v2);
      template_ = temp;
      temp->addEdge(this);
      this->addFacet(facet);
      (v1)->addEdge(this);
      (v2)->addEdge(this);
    }
  }

  // Destructor
  Edge::~Edge()
  {
    mspNode.clear();
    template_ = nullptr;
    mNodes[0] = nullptr;
    mNodes[1] = nullptr;
    facetsConected.clear();
  }

  // Get initial length of the edge
  double Edge::getDist() { return InitialDist; }

  // Get index of the nodes in the edge
  std::pair<uint, uint> Edge::getPairNodes()
  {
    return {mNodes[0]->getIndex(), mNodes[1]->getIndex()};
  }

  // Add facet
  void Edge::addFacet(Facet *facet) { facetsConected.insert(facet); }

  // Erase facet
  void Edge::eraseFacet(Facet *facet) { facetsConected.erase(facet); }

  // Get ptr to the neighbour facets splited by the edge.
  std::set<Facet *> Edge::getFacets() { return facetsConected; }

  // Get ptr to the nodes in the edge.
  std::set<Node *> Edge::getNodes() { return mspNode; }

  // Check if the are no connected facets.
  bool Edge::noFacets() { return (facetsConected.size() == 0); }

  // Delete edge
  void Edge::setBadFlag()
  {
    template_->EraseEdge(this);
    for (std::set<Node *>::iterator it = mspNode.begin(); it != mspNode.end();
         it++)
    {
      (*it)->eraseEdge(this);
    }
  }

  //Comparing functions
  // compare this edge against another edge of the same class
  bool Edge::operator==(const Edge &OtherEdge)
  {
    bool b1 = (this->mNodes[0] == OtherEdge.mNodes[0]) or
              (this->mNodes[0] == OtherEdge.mNodes[1]);
    bool b2 = (this->mNodes[1] == OtherEdge.mNodes[0]) or
              (this->mNodes[1] == OtherEdge.mNodes[1]);
    return b1 && b2;
  }

  // See if this edge links the nodes v1 and v2.
  bool Edge::isEqual(const Node *v1, const Node *v2)
  {
    bool b1 = (this->mNodes[0] == v1) or (this->mNodes[0] == v2);
    bool b2 = (this->mNodes[1] == v1) or (this->mNodes[1] == v2);
    return b1 && b2;
  }
} // namespace defSLAM
