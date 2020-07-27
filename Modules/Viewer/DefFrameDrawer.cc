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

#include "DefFrameDrawer.h"
#include "Tracking.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <mutex>
#include <opencv2/core/eigen.hpp>

namespace defSLAM
{
  using namespace cv;

  /*********
   *  Function to draw a line with 2 colors in opencv
   * inspired in:
   * https://www.rojtberg.net/929/how-to-draw-a-line-interpolating-2-colors-with-opencv/
   * ********/
  void line2(Mat &img, const Point &start, const Point &end, const Scalar &c1,
             const Scalar &c2)
  {
    LineIterator iter(img, start, end, 8);

    for (int i = 0; i < iter.count; i++, iter++)
    {
      double alpha = double(i) / iter.count;
      // note: using img.at<T>(iter.pos()) is faster, but
      // then you have to deal with mat type and channel number yourself
      img(Rect(iter.pos(), Size(1, 1))) = c1 * (1.0 - alpha) + c2 * alpha;
    }
  }

  // Constructor.
  DefFrameDrawer::DefFrameDrawer(Map *pMap) : FrameDrawer(pMap) {}

  // Draw a frame to display. It draw the map, with keypoints and the template if selected.
  // It use the function DrawFrame from FrameDrawer of ORBSLAM_2.
  cv::Mat DefFrameDrawer::DrawFrame()
  {

    cv::Mat im = FrameDrawer::DrawFrame();

    // Set flag to draw template. TODO: It still gives some failures.
    if (false)
    {
      this->DrawTemplate(im);
    }

    return im;
  }

  // Draw the projection of the edges of the template.
  // It is very ilustrative but quite dizzy
  void DefFrameDrawer::DrawTemplate(cv::Mat &im)
  {
    unique_lock<mutex> lock(mMutex);
    std::unique_lock<std::mutex> M(
        static_cast<DefMap *>(mpMap)->MutexUpdating);

    if ((mTcw.rows == 4) and (mTcw.cols == 4) and (mK.rows == 3) and
        (mK.cols == 3))
    {

      std::set<Node *> Nodes =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getNodes();

      vector<cv::Point3f> Nodes3d;
      vector<cv::Point2f> Nodes2d;
      int i(0);
      std::map<uint, uint> M;
      std::vector<Node *> VecNodes;
      for (std::set<Node *>::iterator it = Nodes.begin(); it != Nodes.end();
           it++)
      {
        double x, y, z;
        (*it)->getXYZ(x, y, z);
        Nodes3d.push_back(cv::Point3f(x, y, z));
        M[(*it)->getIndex()] = i;
        VecNodes.push_back(*it);
        i++;
      }

      cv::Mat RotVect;
      cv::Mat rotMat = mTcw(cv::Range(0, 3), cv::Range(0, 3));
      cv::Mat transMat = mTcw(cv::Range(0, 3), cv::Range(3, 4));

      cv::Rodrigues(rotMat, RotVect);
      cv::projectPoints(Nodes3d, RotVect, transMat, mK, cv::noArray(), Nodes2d);
      std::set<Edge *> Edges =
          static_cast<DefMap *>(mpMap)->GetTemplate()->getEdges();
      for (std::set<Edge *>::iterator ite = Edges.begin(); ite != Edges.end();
           ite++)
      {
        std::pair<uint, uint> mo = (*ite)->getPairNodes();
        cv::Scalar S1, S2;
        uint a = M[mo.first];
        uint b = M[mo.second];
        Node *n1 = VecNodes[a];
        Node *n2 = VecNodes[b];

        if (n1->isBoundary())
        {
          S1 = cv::Scalar(1 * 254, 0 * 254, 0.1 * 254);
        }
        else
        {
          if (n1->isViewed())
          {
            S1 = cv::Scalar(0.1 * 254, 0.9 * 254, 0.1 * 254);
          }
          else if (n1->isLocal())
          {
            S1 = cv::Scalar(0.1 * 254, 0.8 * 254, 1.0 * 254);
          }
          else
          {
            S1 = cv::Scalar(0.05 * 254, 0.05 * 254, 0.05 * 254);
          }
        }
        if (n2->isBoundary())
        {
          S2 = cv::Scalar(1 * 254, 0 * 254, 0.1 * 254);
        }
        else
        {
          if (n2->isViewed())
          {
            S2 = cv::Scalar(0.1 * 254, 0.9 * 254, 0.1 * 254);
          }
          else if (n2->isLocal())
          {
            S2 = cv::Scalar(0.1 * 254, 0.8 * 254, 1.0 * 254);
          }
          else
          {
            S2 = cv::Scalar(0.05 * 254, 0.05 * 254, 0.05 * 254);
          }
        }

        try
        {
          if (((Nodes2d[a].x > 0) and (Nodes2d[a].y > 0) and
               (Nodes2d[a].x < im.cols) and (Nodes2d[a].y < im.rows)) or
              ((Nodes2d[b].x > 0) and (Nodes2d[b].y > 0) and
               (Nodes2d[b].x < im.cols) and (Nodes2d[b].y < im.rows)))
            line2(im, Nodes2d[a], Nodes2d[b], S1, S2);
          //         cv::line(im,Nodes2d[a],Nodes2d[b],cv::Scalar(15,15,155),1);
        }
        catch (...)
        {
          std::cout << "Failure drawing template" << std::endl;
          exit(1);
        }
      }
    }
  }

  // Update frame drawer from tracker.
  void DefFrameDrawer::Update(Tracking *pTracker)
  {
    unique_lock<mutex> lock(mMutex);
    pTracker->mImRGB.copyTo(mIm);
    auto points = mpMap->GetReferenceMapPoints();
    mvCurrentLocalMap.clear();
    for (auto pMP : points)
    {
      if (pMP)
      {
        if (pMP->isBad())
          continue;
        if (static_cast<DefMapPoint *>(pMP)->getFacet())
          if (pTracker->mCurrentFrame->isInFrustum(pMP, 0.5))
          {
            cv::KeyPoint kp =
                pTracker->mCurrentFrame->ProjectPoints(pMP->GetWorldPos());
            mvCurrentLocalMap.push_back(std::move(kp));
          }
      }
    }

    mvCurrentKeys = pTracker->mCurrentFrame->mvKeys;
    this->mvCurrentKeysCorr = pTracker->mCurrentFrame->mvKeysUnCorr;
    N = mvCurrentKeys.size();
    this->N2 = mvCurrentKeysCorr.size();

    mvbVO = vector<bool>(N, false);
    mvbMap = vector<bool>(N, false);
    this->mvbMapcorr = vector<bool>(N2, false);

    if (pTracker->mLastProcessedState == Tracking::NOT_INITIALIZED)
    {
      mvIniKeys = pTracker->mInitialFrame.mvKeys;
      mvIniMatches = pTracker->mvIniMatches;
    }
    else if (pTracker->mLastProcessedState == Tracking::OK)
    {
      for (int i = 0; i < N; i++)
      {
        MapPoint *pMP = pTracker->mCurrentFrame->mvpMapPoints[i];
        if (pMP)
        {
          if (!pTracker->mCurrentFrame->mvbOutlier[i])
          {
            if (static_cast<DefMapPoint *>(pMP)->getFacet())
            {
              if (pMP->Observations() > 0)
                mvbMap[i] = true;
              else
                mvbVO[i] = true;
            }
          }
        }
      }
    }
    try
    {

      pTracker->mCurrentFrame->mK.copyTo(mK);
      pTracker->mCurrentFrame->mTcw.copyTo(mTcw);
    }
    catch (...)
    {
      std::cout << "Failure updating position in the frame drawer" << std::endl;
      exit(1);
    }
    mState = static_cast<int>(pTracker->mLastProcessedState);
  }

} // namespace defSLAM
