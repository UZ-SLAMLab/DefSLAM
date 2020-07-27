/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University
*of Zaragoza) && Shaifali Parashar, Adrien Bartoli (Université Clermont Auvergne)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/

#include "MapDrawer.h"
#include "GroundTruthFrame.h"
#include "GroundTruthKeyFrame.h"
#include "KeyFrame.h"
#include "MapPoint.h"
#include <mutex>
#include <pangolin/pangolin.h>
namespace ORB_SLAM2
{
  using defSLAM::GroundTruthFrame;

  MapDrawer::MapDrawer(Map *pMap, const string &strSettingPath) : mpMap(pMap)
  {
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);

    mKeyFrameSize = fSettings["Viewer.KeyFrameSize"];
    mKeyFrameLineWidth = fSettings["Viewer.KeyFrameLineWidth"];
    mGraphLineWidth = fSettings["Viewer.GraphLineWidth"];
    mPointSize = fSettings["Viewer.PointSize"];
    mCameraSize = fSettings["Viewer.CameraSize"];
    mCameraLineWidth = fSettings["Viewer.CameraLineWidth"];
  }

  void MapDrawer::reset() {}

  void MapDrawer::DrawMapPoints()
  {
    const vector<MapPoint *> &vpMPs = mpMap->GetAllMapPoints();
    const vector<MapPoint *> &vpRefMPs = mpMap->GetReferenceMapPoints();

    set<MapPoint *> spRefMPs(vpRefMPs.begin(), vpRefMPs.end());

    if (vpMPs.empty())
      return;

    glPointSize(mPointSize);
    glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);

    for (size_t i = 0, iend = vpMPs.size(); i < iend; i++)
    {
      if (vpMPs[i]->isBad() || spRefMPs.count(vpMPs[i]))
        continue;
      cv::Mat pos = vpMPs[i]->GetWorldPos();
      glVertex3f(pos.at<float>(0), pos.at<float>(1), pos.at<float>(2));
    }
    glEnd();

    glPointSize(mPointSize);

    glBegin(GL_POINTS);

    for (set<MapPoint *>::iterator sit = spRefMPs.begin(), send = spRefMPs.end();
         sit != send; sit++)
    {
      glPointSize(mPointSize);

      if ((*sit)->isBad())
        continue;
      glColor3f(1.0, 0.0, 0.0); // 1.0,0.0,0.0
      cv::Mat pos = (*sit)->GetWorldPos();
      glVertex3f(pos.at<float>(0), pos.at<float>(1), pos.at<float>(2));
    }

    glEnd();
  }

  void MapDrawer::DrawPoints()
  {
    std::unique_lock<std::mutex> K(mPoints);
    glBegin(GL_POINTS);

    glPointSize(5);
    glColor3f(1.0, 0.0, 0.0); // RGB
    for (size_t i = 0, iend = PointsMap.size(); i < iend; i = i + 3)
    {
      //  PointsMap[i+2] << std::endl;
      glVertex3f(PointsMap[i], PointsMap[i + 1],
                 PointsMap[i + 2]); // Red points: All map points non Scaled
    }
    glEnd();

    glBegin(GL_POINTS);

    glPointSize(10);
    glColor3f(0.0, 0.2, 1.0); // RGB
    for (size_t i = 0, iend = PointsRef.size(); i < iend; i = i + 3)
    {
      //  PointsMap[i+2] << std::endl;
      glVertex3f(PointsRef[i], PointsRef[i + 1],
                 PointsRef[i + 2]); // Red points: All map points non Scaled
    }
    glEnd();
  }
  void MapDrawer::UpdatePoints(Frame *pFrame)
  {
    std::unique_lock<std::mutex> K(mPoints);
    // unique_lock<mutex> lock(mMutex);

    // PointsSeen.clear();
    PointsStereoFrame.clear();
    PointsLocal.clear();
    PointsMap.clear();
    if (pFrame->mTcw.empty())
      return;
    std::vector<MapPoint *> MapPoints = pFrame->mvpMapPoints;
    cv::Mat Rcw = pFrame->mTcw.rowRange(0, 3).colRange(0, 3).clone();
    cv::Mat tcw = pFrame->mTcw.rowRange(0, 3).col(3).clone();

    for (size_t i = 0, iend = MapPoints.size(); i < iend; i++)
    {
      if (MapPoints[i])
      {

        float z = pFrame->mvDepth[i];
        if ((z > 0) && (z < 1000)) // 0.7
        {
          /*cv::Mat mono = Rcw*(MapPoints[i]->GetWorldPos())+tcw;
                cv::Mat stereo = (Rcw*(pFrame->UnprojectStereo(i))+tcw);*/
          // Stereo
          cv::Mat mono = MapPoints[i]->GetWorldPos();
          cv::Mat stereo = pFrame->UnprojectStereo(i);

          //                PointsSeen.push_back(mono.at<float>(0));
          //                PointsSeen.push_back(mono.at<float>(1));
          //                PointsSeen.push_back(mono.at<float>(2));

          PointsStereoFrame.push_back(stereo.at<float>(0));
          PointsStereoFrame.push_back(stereo.at<float>(1));
          PointsStereoFrame.push_back(stereo.at<float>(2));
        }
      }
    }

    const vector<MapPoint *> &vpRefMPs = mpMap->GetReferenceMapPoints();
    const vector<MapPoint *> &vpMPs = mpMap->GetAllMapPoints();

    // const vector<MapPoint*> &vpVisibleMPs = mpMap->GetVisibleMapPoints();
    set<MapPoint *> spVisibleMPs(vpMPs.begin(), vpMPs.end());

    set<MapPoint *> spRefMPs(vpRefMPs.begin(), vpRefMPs.end());

    for (set<MapPoint *>::iterator sit = spRefMPs.begin(), send = spRefMPs.end();
         sit != send; sit++)
    {
      if ((*sit)->isBad())
        continue;
      cv::Mat pos = (*sit)->GetWorldPos();
      PointsLocal.push_back(pos.at<float>(0));
      PointsLocal.push_back(pos.at<float>(1));
      PointsLocal.push_back(pos.at<float>(2));
    }

    for (set<MapPoint *>::iterator sit = spVisibleMPs.begin(),
                                   send = spVisibleMPs.end();
         sit != send; sit++)
    {
      if ((*sit)->isBad())
        continue;
      cv::Mat pos = (*sit)->GetWorldPos();
      PointsMap.push_back(pos.at<float>(0));
      PointsMap.push_back(pos.at<float>(1));
      PointsMap.push_back(pos.at<float>(2));
      // pointsPlaying++;
    }
    //    cout << "Visible Points: " << spVisibleMPs.size() << endl;
    //    cout << "Reference points: " << vpRefMPs.size() << endl;
    //    cout << "AllMapPoints: " << vpMPs.size() << endl;
    //    cout << "MapPoints in facet: " << pointsPlaying << endl;
  }

  void MapDrawer::UpdatePoints(Frame *pFrame, float s)
  {
    std::unique_lock<std::mutex> K(mPoints);
    // unique_lock<mutex> lock(mMutex);

    PointsSeen.clear();
    PointsMono.clear();
    PointsGT.clear();
    PointsMap.clear();

    std::vector<MapPoint *> AllMapPoints = mpMap->GetAllMapPoints();

    std::vector<std::vector<float>> LocalPoints =
        static_cast<GroundTruthFrame *>(pFrame)->getPosMono();
    std::vector<std::vector<float>> StereoPoints =
        static_cast<GroundTruthFrame *>(pFrame)->getPosStereo();

    cv::Mat Rwc = pFrame->GetRotationInverse();
    cv::Mat twc = pFrame->GetCameraCenter();

    for (size_t i = 0, iend = LocalPoints.size(); i < iend; i++)
    {
      cv::Mat mono_c = (cv::Mat_<float>(3, 1) << LocalPoints[i][0],
                        LocalPoints[i][1], LocalPoints[i][2]);
      cv::Mat stereo_c = (cv::Mat_<float>(3, 1) << StereoPoints[i][0],
                          StereoPoints[i][1], StereoPoints[i][2]);

      cv::Mat mono_c_scaled = mono_c;

      cv::Mat mono_w = Rwc * mono_c + twc;
      cv::Mat mono_w_scaled = Rwc * mono_c_scaled + twc;
      cv::Mat stereo_w = Rwc * stereo_c / s + twc;

      PointsMono.push_back(mono_w.at<float>(0));
      PointsMono.push_back(mono_w.at<float>(1));
      PointsMono.push_back(mono_w.at<float>(2));

      PointsSeen.push_back(mono_w_scaled.at<float>(0));
      PointsSeen.push_back(mono_w_scaled.at<float>(1));
      PointsSeen.push_back(mono_w_scaled.at<float>(2));

      PointsGT.push_back(stereo_w.at<float>(0));
      PointsGT.push_back(stereo_w.at<float>(1));
      PointsGT.push_back(stereo_w.at<float>(2));
    }
    for (uint i(0); i < AllMapPoints.size(); i++)
    {
      MapPoint *pMP = AllMapPoints[i];
      if (!pMP)
        continue;
      if (pMP->isBad())
        continue;
      cv::Mat posMP = pMP->GetWorldPos();
      PointsMap.push_back(posMP.at<float>(0));
      PointsMap.push_back(posMP.at<float>(1));
      PointsMap.push_back(posMP.at<float>(2));
    }
  }

  void MapDrawer::DrawGTPoints()
  {
    std::unique_lock<std::mutex> K(mPoints);

    for (size_t i = 0, iend = PointsSeen.size(); i < iend; i = i + 3)
    {
      glPointSize(5);
      glBegin(GL_POINTS);

      // Sometimes stereo is an empty cv::Mat-->need to be solved
      if (PointsGT[i + 2] > 0.0)
      {
        double red[3] = {62.0 / 255.0, 37.0 / 255, 41.0 / 255};
        glColor3f(red[0], red[1], red[2]); // RGB
        glVertex3f(PointsMono[i], PointsMono[i + 1], PointsMono[i + 2]);
      }

      glEnd();
    }

    for (size_t i = 0, iend = PointsSeen.size(); i < iend; i = i + 3)
    {
      glVertex3f(PointsSeen[i], PointsSeen[i + 1],
                 PointsSeen[i + 2]); // Green points: mono points scaled
    }
    glEnd();

    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0); // BGR

    glEnd();

    glPointSize(3);

    glBegin(GL_POINTS);
    glColor3f(1.0, 1.0, 0.0);

    for (size_t i = 0, iend = PointsGT.size(); i < iend; i = i + 3)
    {
      glVertex3f(PointsGT[i], PointsGT[i + 1],
                 PointsGT[i + 2]); // Yellow points: stereo points (GT)
    }
    glEnd();

    // Draw error lines
    float lambda = 1.0;
    for (size_t i = 0, iend = PointsSeen.size(); i < iend; i = i + 3)
    {
      glLineWidth(mKeyFrameLineWidth); //++3
      glColor3f(0.0, 0.0, 1.0);
      glBegin(GL_LINES);

      // Sometimes stereo is an empty cv::Mat-->need to be solved
      if (PointsGT[i + 2] > 0.0)
      {
        glVertex3f(PointsSeen[i], PointsSeen[i + 1], PointsSeen[i + 2]);
        float v1, v2, v3;
        v1 = PointsSeen[i] + lambda * (PointsGT[i] - PointsSeen[i]);
        v2 = PointsSeen[i + 1] + lambda * (PointsGT[i + 1] - PointsSeen[i + 1]);
        v3 = PointsSeen[i + 2] + lambda * (PointsGT[i + 2] - PointsSeen[i + 2]);
        glVertex3f(v1, v2, v3);
      }

      glEnd();
    }
  }

  void MapDrawer::DrawError()
  {
    const vector<MapPoint *> &vpMPs = mpMap->GetAllMapPoints();
    const vector<MapPoint *> &vpRefMPs = mpMap->GetReferenceMapPoints();

    set<MapPoint *> spRefMPs(vpRefMPs.begin(), vpRefMPs.end());

    if (vpMPs.empty())
      return;

    for (size_t i = 0, iend = vpMPs.size(); i < iend; i++)
    {
      if (vpMPs[i]->isBad() || spRefMPs.count(vpMPs[i]))
        continue;
      cv::Mat pos = vpMPs[i]->GetWorldPos();
      glBegin(GL_POINT);
      glVertex3f(pos.at<float>(0), pos.at<float>(1), pos.at<float>(2));
      glEnd();
    }
  }

  void MapDrawer::DrawKeyFrames(const bool bDrawKF, const bool bDrawGraph)
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
        glColor3f(0.0f, 0.0f, 1.0f);
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

  void MapDrawer::DrawCurrentCamera(pangolin::OpenGlMatrix &Twc)
  {
    const float &w = mCameraSize;
    const float h = w * 0.75;
    const float z = w * 0.6;

    glPushMatrix();

#ifdef HAVE_GLES
    glMultMatrixf(Twc.m);
#else
    glMultMatrixd(Twc.m);
#endif

    glLineWidth(mCameraLineWidth);
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

  void MapDrawer::SetCurrentCameraPose(const cv::Mat &Tcw)
  {
    unique_lock<mutex> lock(mMutexCamera);
    mCameraPose = Tcw.clone();
  }

  void MapDrawer::GetCurrentOpenGLCameraMatrix(pangolin::OpenGlMatrix &M)
  {
    if (!mCameraPose.empty())
    {
      cv::Mat Rwc(3, 3, CV_32F);
      cv::Mat twc(3, 1, CV_32F);
      {
        unique_lock<mutex> lock(mMutexCamera);
        Rwc = mCameraPose.rowRange(0, 3).colRange(0, 3).t();
        twc = -Rwc * mCameraPose.rowRange(0, 3).col(3);
      }

      M.m[0] = Rwc.at<float>(0, 0);
      M.m[1] = Rwc.at<float>(1, 0);
      M.m[2] = Rwc.at<float>(2, 0);
      M.m[3] = 0.0;

      M.m[4] = Rwc.at<float>(0, 1);
      M.m[5] = Rwc.at<float>(1, 1);
      M.m[6] = Rwc.at<float>(2, 1);
      M.m[7] = 0.0;

      M.m[8] = Rwc.at<float>(0, 2);
      M.m[9] = Rwc.at<float>(1, 2);
      M.m[10] = Rwc.at<float>(2, 2);
      M.m[11] = 0.0;

      M.m[12] = twc.at<float>(0);
      M.m[13] = twc.at<float>(1);
      M.m[14] = twc.at<float>(2);
      M.m[15] = 1.0;
    }
    else
      M.SetIdentity();
  }

} // namespace ORB_SLAM2
