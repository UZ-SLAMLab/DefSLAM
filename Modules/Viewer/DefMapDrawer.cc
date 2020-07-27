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

#include "DefMapDrawer.h"
#include "DefKeyFrame.h"
#include "DefMap.h"
#include "DefMapPoint.h"
#include "Edge.h"
#include "Facet.h"
#include "GroundTruthFrame.h"
#include "GroundTruthKeyFrame.h"
#include "KeyFrame.h"
#include "MeshDrawer.h"
#include "Node.h"
#include "Template.h"
#include <mutex>
#include <pangolin/pangolin.h>

namespace defSLAM
{
  // Constructor.
  DefMapDrawer::DefMapDrawer(Map *pMap,
                             const string &strSettingPath)
      : MapDrawer(pMap, strSettingPath), MeshDrawers(nullptr)
  {
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
  }

  // Clear template history.
  void DefMapDrawer::reset() { this->MeshDrawershist.clear(); }

  // Updates from the tracking.
  void DefMapDrawer::updateTemplate()
  {
    std::unique_lock<std::mutex> M(this->mTemplate);
    if (static_cast<DefMap *>(mpMap)->GetTemplate())
    {
      std::unique_lock<std::mutex> lc(
          static_cast<DefMap *>(mpMap)->GetTemplate()->mutexDrawer_);

      if (this->MeshDrawers)
      {
        delete this->MeshDrawers;
      }

      this->MeshDrawers = new MeshDrawer();
      cv::Mat imRGB =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getTexture();

      if ((Texture) && (!imRGB.empty()))
      {
        this->MeshDrawers->addTextureImage(imRGB);
      }
      else
      {
        cv::Mat imRGBNoTex(1024, 1024, CV_8UC3, cv::Scalar(155, 155, 155));
        this->MeshDrawers->addTextureImage(imRGBNoTex);
      }

      std::set<Node *> Nodes =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getNodes();
      std::vector<Node *> NodesVec(Nodes.begin(), Nodes.end());

      for (uint i(0); i < NodesVec.size(); i++)
      {
        std::vector<double> pos;
        std::vector<double> proj;
        pos.push_back(NodesVec[i]->x);
        pos.push_back(NodesVec[i]->y);
        pos.push_back(NodesVec[i]->z);
        proj.push_back(NodesVec[i]->proju / imRGB.cols);
        proj.push_back(NodesVec[i]->projv / imRGB.rows);
        NodesVec[i]->posVecDrawer = i;

        if (static_cast<Node *>(NodesVec[i])->isBoundary())
        {
          this->MeshDrawers->addNode(pos, proj, 0);
        }
        else
        {
          if (static_cast<Node *>(NodesVec[i])->isViewed())
          {
            this->MeshDrawers->addNode(pos, proj, 1);
          }
          else if (static_cast<Node *>(NodesVec[i])->isLocal())
          {
            this->MeshDrawers->addNode(pos, proj, 2);
          }
          else
          {
            this->MeshDrawers->addNode(pos, proj, 4);
          }
        }
      }
      std::set<Facet *> Facets =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getFacets();

      for (std::set<Facet *>::iterator itf = Facets.begin(); itf != Facets.end();
           itf++)
      {
        std::vector<int> fac;
        std::set<Node *> Nodes = (*itf)->getNodes();
        for (std::set<Node *>::iterator it = Nodes.begin(); it != Nodes.end();
             it++)
        {
          fac.push_back((*it)->posVecDrawer);
        }

        this->MeshDrawers->addFacet(fac);
      }

      // Draw Edges
      std::set<Edge *> Edges =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getEdges();
      for (std::set<Edge *>::iterator ite = Edges.begin(); ite != Edges.end();
           ite++)
      {
        std::set<Node *> Nodesedge = (*ite)->getNodes();
        std::vector<int> edge;
        for (std::set<Node *>::iterator it = Nodesedge.begin();
             it != Nodesedge.end(); it++)
        {
          edge.push_back((*it)->posVecDrawer);
        }
        this->MeshDrawers->addEdge(edge);
      }
    }
  }

  // Update template at rest.
  void DefMapDrawer::updateTemplateAtRest()
  {
    std::unique_lock<std::mutex> M(this->mTemplatehist);
    if (static_cast<DefMap *>(mpMap)->ThereIsATemplate)
    {
      auto defMap = static_cast<DefMap *>(mpMap);
      auto templateToDraw = defMap->GetTemplate();
      auto refKf = templateToDraw->kf;

      std::unique_lock<std::mutex> lc(templateToDraw->mutexDrawer_);
      this->MeshDrawershist[refKf] = unique_ptr<MeshDrawer>(new MeshDrawer());
      MeshDrawer *temp = this->MeshDrawershist[refKf].get();
      cv::Mat imRGB = templateToDraw->getTexture();

      if ((Texture) && (!imRGB.empty()))
      {
        temp->addTextureImage(imRGB);
      }
      else
      {
        cv::Mat imRGBNoTex(1024, 1024, CV_8UC3, cv::Scalar(155, 155, 155));
        temp->addTextureImage(imRGBNoTex);
      }

      std::set<Node *> Nodes = templateToDraw->getNodes();
      std::vector<Node *> NodesVec(Nodes.begin(), Nodes.end());

      for (uint i(0); i < NodesVec.size(); i++)
      {
        std::vector<double> pos;
        std::vector<double> proj;
        pos.push_back(NodesVec[i]->xO);
        pos.push_back(NodesVec[i]->yO);
        pos.push_back(NodesVec[i]->zO);
        proj.push_back(NodesVec[i]->proju / imRGB.cols);
        proj.push_back(NodesVec[i]->projv / imRGB.rows);
        NodesVec[i]->posVecDrawer = i;
        temp->addNode(pos, proj, 4);
      }
      std::set<Facet *> Facets = templateToDraw->getFacets();

      for (std::set<Facet *>::iterator itf = Facets.begin(); itf != Facets.end();
           itf++)
      {
        std::vector<int> fac;
        std::set<Node *> Nodes = (*itf)->getNodes();
        for (std::set<Node *>::iterator it = Nodes.begin(); it != Nodes.end();
             it++)
        {
          fac.push_back((*it)->posVecDrawer);
        }
        temp->addFacet(fac);
      }

      // Draw Edges
      std::set<Edge *> Edges = defMap->GetTemplate()->getEdges();
      for (std::set<Edge *>::iterator ite = Edges.begin(); ite != Edges.end();
           ite++)
      {
        std::set<Node *> Nodesedge = (*ite)->getNodes();
        std::vector<int> edge;
        for (std::set<Node *>::iterator it = Nodesedge.begin();
             it != Nodesedge.end(); it++)
        {
          edge.push_back((*it)->posVecDrawer);
        }
        temp->addEdge(edge);
      }
      PointsRef.clear();
      for (uint i(0); i < refKf->mvKeys.size(); i++)
      {
        MapPoint *pMP = refKf->GetMapPoint(i);
        if (pMP)
        {
          if (pMP->isBad())
            continue;

          if (static_cast<DefMapPoint *>(pMP)->getFacet())
          {
            if (!pMP->covNorm)
              continue;
            cv::Mat x3Dy = static_cast<DefMapPoint *>(pMP)
                               ->PosesKeyframes[refKf]
                               .clone();
            if (x3Dy.empty())
              continue;
            PointsRef.push_back(x3Dy.at<float>(0, 0));
            PointsRef.push_back(x3Dy.at<float>(1, 0));
            PointsRef.push_back(x3Dy.at<float>(2, 0));
          }
        }
      }
    }
  }

  // Draw the template through meshdrawer.
  void DefMapDrawer::DrawTemplate()
  {
    std::unique_lock<std::mutex> M(this->mTemplate);
    if (this->MeshDrawers)
    {
      this->MeshDrawers->drawMesh();
    }
  }

  // Draw the shape-at-rest of the template.
  void DefMapDrawer::DrawTemplateAtRest(uint o)
  {
    std::unique_lock<std::mutex> M(
        static_cast<DefMap *>(mpMap)->MutexUpdating);

    if (static_cast<DefMap *>(mpMap)->ThereIsATemplate)
    {
      std::unique_lock<std::mutex> lc(
          static_cast<DefMap *>(mpMap)->GetTemplate()->mutexDrawer_);

      std::set<Facet *> Facets =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getFacets();

      if (Texture)
      {
        glEnable(GL_TEXTURE_2D);

        for (std::map<KeyFrame *, std::set<Facet *>>::iterator itKf =
                 Facet::TexturesDataset.begin();
             itKf != Facet::TexturesDataset.end(); itKf++)
        {
          GLuint texID;
          glGenTextures(1, &texID);
          glBindTexture(GL_TEXTURE_2D, texID);
          glPixelStoref(GL_UNPACK_ALIGNMENT, 1);
          glLineWidth(1);

          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

          glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, itKf->first->RGBimage.cols,
                       itKf->first->RGBimage.rows, 0, GL_BGR, GL_UNSIGNED_BYTE,
                       itKf->first->RGBimage.ptr());

          for (std::set<Facet *>::iterator itf = itKf->second.begin();
               itf != itKf->second.end(); itf++)
          {
            glBegin(GL_TRIANGLES);
            glColor3f(1, 1, 1); // set global color to white, otherwise this color
                                // will be (somehow) added to the texture

            // glColor3d((*itf)->R,(*itf)->G,(*itf)->B);
            std::set<Node *> Nodes = (*itf)->getNodes();
            for (std::set<Node *>::iterator it = Nodes.begin(); it != Nodes.end();
                 it++)
            {
              cv::Point UV = (*itf)->texture.second[(*it)];
              glTexCoord2f(((float)UV.x) / itKf->first->RGBimage.cols,
                           ((float)UV.y) / itKf->first->RGBimage.rows);
              glVertex3f((*it)->xO, (*it)->yO, (*it)->zO);
            }
          }
          glEnd();
          glDeleteTextures(1, &texID); // TODO load all textures once and only
                                       // remove when it change
        }

        glDisable(GL_TEXTURE_2D);

        for (std::set<Facet *>::iterator itf = Facets.begin();
             itf != Facets.end(); itf++)
        {
          // Facets
          if (!(*itf)->texture.first)
          {
            /*  glBegin(GL_TRIANGLES);
            glColor3d((*itf)->R,(*itf)->G,(*itf)->B);
            std::array<Node*, 3>  Nodes((*itf)->getNodesArray());

            for(uint i(0);i<3;i++){
                glVertex3f((Nodes[i])->x,(Nodes[i])->y,(Nodes[i])->z);
                j++;
            }*/
          }
        }
      }

      // Draw Edges
      std::set<Edge *> Edges =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getEdges();
      for (std::set<Edge *>::iterator ite = Edges.begin(); ite != Edges.end();
           ite++)
      {
        std::set<Node *> Nodesedge = (*ite)->getNodes();
        glBegin(GL_LINES);

        for (std::set<Node *>::iterator it = Nodesedge.begin();
             it != Nodesedge.end(); it++)
        {
          if (static_cast<Node *>(*it)->isBoundary())
          {
            glColor3f(0.1, 0, 1);
            glVertex3f((*it)->xO, (*it)->yO, (*it)->zO);
          }
          else
          {
            if (static_cast<Node *>(*it)->isViewed())
            {
              glColor3f(0.1, 0.9, 0.1);
              glVertex3f((*it)->xO, (*it)->yO, (*it)->zO);
            }
            else if (static_cast<Node *>(*it)->isLocal())
            {
              glColor3f(1.0, 0.8, 0.1);
              glVertex3f((*it)->xO, (*it)->yO, (*it)->zO);
            }
            else
            {
              glColor3f(0.05, 0.05, 0.05);
              glVertex3f((*it)->xO, (*it)->yO, (*it)->zO);
            }
          }
        }
        glEnd();
      }
    }
  }

  // Draw the templates through the execution.
  void DefMapDrawer::DrawTemplatehist()
  {
    std::unique_lock<std::mutex> M(this->mTemplatehist);
    for (auto &it : this->MeshDrawershist)
    {
      double alpha(0.99);
      it.second->drawMesh(alpha, false);
    }
  }

  // Draw template with texture.
  void DefMapDrawer::ShowTexture(bool C)
  {
    std::unique_lock<std::mutex> mL(mTexture);
    Texture = C;
  }

  // Draw the keyframes. Blue for ordinary keyframes and red for
  // anchor keyframe.
  void DefMapDrawer::DrawKeyFrames(const bool bDrawKF,
                                   const bool bDrawGraph)
  {
    const float &w = mKeyFrameSize;
    const float h = w * 0.75;
    const float z = w * 0.6;

    const vector<KeyFrame *> vpKFs = mpMap->GetAllKeyFrames();

    if (bDrawKF)
    {
      for (size_t i = 0; i < vpKFs.size(); i++)
      {
        KeyFrame *pKF = vpKFs[i];
        cv::Mat Twc = pKF->GetPoseInverse().t();

        glPushMatrix();

        glMultMatrixf(Twc.ptr<GLfloat>(0));

        glLineWidth(mKeyFrameLineWidth);
        if (static_cast<DefKeyFrame *>(pKF)->templateAssigned())
          glColor3f(0.0f, 0.0f, 1.0f);
        else
          glColor3f(1.0f, 0.0f, 0.0f);

        if ((pKF) == static_cast<DefMap *>(mpMap)->getLastKFTemplate())
          glColor3f(0.0f, 1.0f, 0.0f);

        glBegin(GL_LINES);
        glVertex3f(0, 0, 0);
        glVertex3f(w, h, z);
        glVertex3f(0, 0, 0);
        glVertex3f(w, -h, z);
        glVertex3f(0, 0, 0);
        glVertex3f(-w, -h, z);
        glVertex3f(0, 0, 0);
        glVertex3f(-w, h, z);

        glVertex3f(w, h, z);
        glVertex3f(w, -h, z);

        glVertex3f(-w, h, z);
        glVertex3f(-w, -h, z);

        glVertex3f(-w, h, z);
        glVertex3f(w, h, z);

        glVertex3f(-w, -h, z);
        glVertex3f(w, -h, z);
        glEnd();

        glPopMatrix();
      }
    }

    if (bDrawGraph)
    {
      glLineWidth(mGraphLineWidth);
      glColor4f(0.0f, 1.0f, 0.0f, 0.6f);
      glBegin(GL_LINES);

      for (size_t i = 0; i < vpKFs.size(); i++)
      {
        // Covisibility Graph
        const vector<KeyFrame *> vCovKFs = vpKFs[i]->GetCovisiblesByWeight(100);
        cv::Mat Ow = vpKFs[i]->GetCameraCenter();
        if (!vCovKFs.empty())
        {
          for (vector<KeyFrame *>::const_iterator vit = vCovKFs.begin(),
                                                  vend = vCovKFs.end();
               vit != vend; vit++)
          {
            if ((*vit)->mnId < vpKFs[i]->mnId)
              continue;
            cv::Mat Ow2 = (*vit)->GetCameraCenter();
            glVertex3f(Ow.at<float>(0), Ow.at<float>(1), Ow.at<float>(2));
            glVertex3f(Ow2.at<float>(0), Ow2.at<float>(1), Ow2.at<float>(2));
          }
        }

        // Spanning tree
        KeyFrame *pParent = vpKFs[i]->GetParent();
        if (pParent)
        {
          cv::Mat Owp = pParent->GetCameraCenter();
          glVertex3f(Ow.at<float>(0), Ow.at<float>(1), Ow.at<float>(2));
          glVertex3f(Owp.at<float>(0), Owp.at<float>(1), Owp.at<float>(2));
        }

        // Loops
        set<KeyFrame *> sLoopKFs = vpKFs[i]->GetLoopEdges();
        for (set<KeyFrame *>::iterator sit = sLoopKFs.begin(),
                                       send = sLoopKFs.end();
             sit != send; sit++)
        {
          if ((*sit)->mnId < vpKFs[i]->mnId)
            continue;
          cv::Mat Owl = (*sit)->GetCameraCenter();
          glVertex3f(Ow.at<float>(0), Ow.at<float>(1), Ow.at<float>(2));
          glVertex3f(Owl.at<float>(0), Owl.at<float>(1), Owl.at<float>(2));
        }
      }

      glEnd();
    }
  }
} // namespace defSLAM
