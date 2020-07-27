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

#include <Facet.h>
#include <mutex>

namespace defSLAM
{
  // Texture of the facet in the keyframe created.
  std::map<KeyFrame *, std::set<Facet *>> Facet::TexturesDataset;

  // Constructor. It create a facet from three ptr to nodes for the
  // template temp.
  Facet::Facet(Node *v1, Node *v2, Node *v3, Template *temp) : template_(temp)
  {
    Nodes.insert(v1);
    (v1)->addFacet(this);
    Nodes.insert(v2);
    (v2)->addFacet(this);
    Nodes.insert(v3);
    (v3)->addFacet(this);
    mNodes[0] = v1;
    mNodes[1] = v2;
    mNodes[2] = v3;

    bool repeated[3] = {false, false, false};
    for (auto edg : template_->getEdges())
    {
      repeated[0] = repeated[0] or edg->isEqual(v1, v2);
      repeated[1] = repeated[1] or edg->isEqual(v2, v3);
      repeated[2] = repeated[2] or edg->isEqual(v1, v3);
    }
    if (!repeated[0])
      new Edge(v1, v2, this, temp);
    if (!repeated[1])
      new Edge(v2, v3, this, temp);
    if (!repeated[2])
      new Edge(v1, v3, this, temp);

    texture.first = static_cast<KeyFrame *>(nullptr);
    R = 0.5;
    G = 0.5;
    B = 0.5;
  }

  // Destructor
  Facet::~Facet()
  {
    Nodes.clear();
    mNodes[0] = static_cast<Node *>(nullptr);
    mNodes[1] = static_cast<Node *>(nullptr);
    mNodes[2] = static_cast<Node *>(nullptr);
    edges_.clear();
    mPoints.clear();
    template_ = static_cast<Template *>(nullptr);
  }

  // get nodes in a facet as a set.
  std::set<Node *> Facet::getNodes()
  {
    return Nodes;
  }

  // get nodes in a facet as an array.
  std::array<Node *, 3> Facet::getNodesArray() { return mNodes; }

  // Check there are observations in the facet.
  bool Facet::noObservations() { return (mPoints.size() == 0); }

  // Add map point observed in the facet.
  void Facet::addMapPoint(MapPoint *mp) { mPoints.insert(mp); }

  // Erase map point observed in the facet.
  void Facet::eraseMapPoint(MapPoint *mp) { mPoints.erase(mp); }

  // Get edges of the facet.
  std::set<Edge *> Facet::getEdges() { return edges_; }

  // Remove facet
  void Facet::setBadFlag()
  {
    template_->eraseFacet(this);

    std::set<Node *> MapNodes = template_->getNodes();
    for (std::set<Node *>::iterator it = MapNodes.begin(); it != MapNodes.end();
         it++)
    {
      (*it)->eraseFacet(this);
      if ((*it)->noFacets())
        (*it)->setBadFlag();
    }
  }

  // get map points in the facet.
  std::set<MapPoint *> Facet::getMapPoints() { return mPoints; }

  // Get texture for the image in the keyframe given
  void Facet::getTextureCoordinates(KeyFrame *KF)
  {
    // Loop over the Keyframes to see a proper one
    texture.first = nullptr;
    vector<cv::Point3f> Nodes3d;
    vector<cv::Point2f> Nodes2d;

    for (std::set<Node *>::iterator itn = Nodes.begin(); itn != Nodes.end();
         itn++)
    {
      double x, y, z;
      (*itn)->getXYZ(x, y, z);
      Nodes3d.push_back(cv::Point3f(float(x), float(y), float(z)));
    }
    cv::Mat mTcw = KF->GetPose();

    cv::Mat mK = KF->mK.clone();
    cv::Mat RotVect;
    cv::Mat rotMat = mTcw(cv::Range(0, 3), cv::Range(0, 3));
    cv::Mat transMat = mTcw(cv::Range(0, 3), cv::Range(3, 4));
    cv::Rodrigues(rotMat, RotVect);
    cv::projectPoints(Nodes3d, RotVect, transMat, mK, cv::noArray(), Nodes2d);

    bool IsIn(true);
    uint i(0);
    std::map<Node *, cv::Point> points;
    for (std::set<Node *>::iterator itn = Nodes.begin(); itn != Nodes.end();
         itn++)
    {
      points[*itn] = Nodes2d[i];
      (*itn)->proju = Nodes2d[i].x;
      (*itn)->projv = Nodes2d[i].y;
      if (!((Nodes2d[i].x < 640) && (Nodes2d[i].y < 480) && (Nodes2d[i].x > 0) &&
            (Nodes2d[i].y > 0)))
      {
        IsIn = false;
      }
      i++;
    }

    if (IsIn)
    {
      texture.first = KF;
      TexturesDataset[KF].insert(this);
      texture.second = points;
    }
  }
} // namespace defSLAM
