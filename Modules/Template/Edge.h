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

#ifndef EDGE_H
#define EDGE_H

#include <Template.h>
#include <Node.h>
#include <Facet.h>
#include <set>
#include <mutex>

namespace defSLAM
{

  class Template;
  class Node;
  class Facet;

  class Edge
  {
  public:
    // Constructor. Initialize an edge from two nodes, a facet
    // and the template
    Edge(Node *, Node *, Facet *, Template *);

    // Destructor
    ~Edge();

    // Get initial length of the edge
    double getDist();

    // Get index of the nodes in the edge
    std::pair<uint, uint> getPairNodes();

    // Add facet
    void addFacet(Facet *facet);

    // Erase facet
    void eraseFacet(Facet *facet);

    // get ptr to the neighbour facets splited by the edge.
    std::set<Facet *> getFacets();

    // get ptr to the nodes in the edge.
    std::set<Node *> getNodes();

    // Check if the are no connected facets.
    bool noFacets();

    // Delete this edge
    void setBadFlag();

    // returns template that contains this edge.
    Template *getTemplate() { return template_; }

    //Comparing functions
    // compare this edge against another edge of the same class
    bool operator==(const Edge &);

    // See if this edge links the nodes v1 and v2.
    bool isEqual(const Node *v1, const Node *v2);

  private:
    Template *template_;      // Template that contains this edge.
    std::set<Node *> mspNode; // Set of ptr to the nodes linked by this edge.

    Node *mNodes[2];                  //
    std::set<Facet *> facetsConected; // set of nodes.
    double InitialDist;               // Distance when the edge was initialized, i.e. shape-at-rest.
  };

} // namespace defSLAM

#endif // EDGE_H
