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

#include "Viewer.h"
#include <iomanip>
#include <mutex>
#include <pangolin/pangolin.h>
//#include <unistd.h>
namespace ORB_SLAM2
{

  Viewer::Viewer(System *pSystem, FrameDrawer *pFrameDrawer,
                 MapDrawer *pMapDrawer, Tracking *pTracking,
                 const string &strSettingPath)
      : mpSystem(pSystem), mpFrameDrawer(pFrameDrawer), mpMapDrawer(pMapDrawer),
        mpTracker(pTracking), mbFinishRequested(false), mbFinished(true),
        mbStopped(true), mbStopRequested(false), mbSaveResults(false),
        next(false)
  {
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);

    float fps = fSettings["Camera.fps"];
    if (fps < 1)
      fps = 30;
    mT = 1e3 / 30;

    mImageWidth = fSettings["Camera.width"];
    mImageHeight = fSettings["Camera.height"];
    if (mImageWidth < 1 || mImageHeight < 1)
    {
      mImageWidth = 640;
      mImageHeight = 480;
    }

    mViewpointX = fSettings["Viewer.ViewpointX"];
    mViewpointY = fSettings["Viewer.ViewpointY"];
    mViewpointZ = fSettings["Viewer.ViewpointZ"];
    mViewpointF = fSettings["Viewer.ViewpointF"];
    double SaveResults = fSettings["Viewer.SaveResults"];
    mbSaveResults = (unsigned int)(SaveResults);
  }

  void Viewer::Run()
  {
    mbFinished = false;
    mbStopped = false;

    pangolin::CreateWindowAndBind("ORBSLAM2: Map Viewer", 175 + 640, 480);

    // 3D Mouse handler requires depth testing to be enabled
    glEnable(GL_DEPTH_TEST);

    // Issue specific OpenGl we might need
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    pangolin::CreatePanel("menu").SetBounds(0.0, 1.0, 0.0,
                                            pangolin::Attach::Pix(175));
    pangolin::Var<bool> menuFollowCamera("menu.Follow Camera", true, true);
    pangolin::Var<bool> menuShowPoints("menu.Show Points", true, true);
    pangolin::Var<bool> menuShowKeyFrames("menu.Show KeyFrames", false, true);
    pangolin::Var<bool> menuShowGraph("menu.Show Graph", false, true);
    pangolin::Var<bool> menuShowError("menu.Show Error", false, true);
    pangolin::Var<bool> menuAutoPlay("menu.Autoplay", true, true);

    // Added for non-rigid method
    pangolin::Var<bool> menuLocalizationMode("menu.Localization Mode", false,
                                             true);

    pangolin::Var<bool> menuNext("menu.Next", false, false);
    pangolin::Var<bool> menuReset("menu.Reset", false, false);

    // Define Camera Render Object (for view / scene browsing)
    pangolin::OpenGlMatrix proj = pangolin::ProjectionMatrix(
        1024, 768, mViewpointF, mViewpointF, 512, 389, 0.1, 1000);
    //  pangolin::OpenGlMatrix proj2 =
    //  pangolin::ProjectionMatrix(1024,768,mViewpointF,mViewpointF,512,389,0.1,1000);

    pangolin::OpenGlRenderState s_cam(
        proj, pangolin::ModelViewLookAt(mViewpointX, mViewpointY, mViewpointZ,
                                        (mViewpointX + (3)), mViewpointY + 0.1,
                                        mViewpointZ + 2, 0.0, -1, 0.0));

    // Add named OpenGL viewport to window and provide 3D Handler
    // Output 3D
    pangolin::View &d_cam1 =
        pangolin::CreateDisplay()
            .SetAspect(640.0f / 480.0f)
            .SetBounds(0.0, 1.0, pangolin::Attach::Pix(175), 1, -1024.0f / 768.0f)
            .SetHandler(new pangolin::Handler3D(s_cam));

    pangolin::OpenGlMatrix Twc;
    Twc.SetIdentity();

    cv::namedWindow("ORBSLAM2: Current Frame");

    bool bFollow = true;
    bool bLocalizationMode = false;

    while (1)
    {
        this_thread::sleep_for(chrono::microseconds((size_t)(mT * 1000)));
      // sleep(1);

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      mpMapDrawer->GetCurrentOpenGLCameraMatrix(Twc);

      if (menuFollowCamera && bFollow)
      {
        s_cam.Follow(Twc);
      }
      else if (menuFollowCamera && !bFollow)
      {
        s_cam.SetModelViewMatrix(pangolin::ModelViewLookAt(
            mViewpointX, mViewpointY, mViewpointZ, 0, 0, 0, 0.0, -1.0, 0.0));
        s_cam.Follow(Twc);
        bFollow = true;
      }
      else if (!menuFollowCamera && bFollow)
      {
        bFollow = false;
      }

      if (menuLocalizationMode && !bLocalizationMode)
      {
        mpSystem->ActivateLocalizationMode();
        bLocalizationMode = true;
      }
      else if (!menuLocalizationMode && bLocalizationMode)
      {
        mpSystem->DeactivateLocalizationMode();
        bLocalizationMode = false;
      }

      d_cam1.Activate(s_cam);

      glClearColor(0.5f, 0.5f, 0.5f, 0.f);
      mpMapDrawer->DrawCurrentCamera(Twc);
      if (menuShowKeyFrames || menuShowGraph)
        mpMapDrawer->DrawKeyFrames(menuShowKeyFrames, menuShowGraph);
      if (menuShowPoints)
        mpMapDrawer->DrawPoints();

      if (menuShowError)
        mpMapDrawer->DrawGTPoints();

      if (mbSaveResults)
      {
        std::ostringstream out;
        out << std::internal << std::setfill('0') << std::setw(5)
            << (unsigned int)(timestamp);
        d_cam1.SaveOnRender("3D" + out.str());
      }
      cv::Mat im = mpFrameDrawer->DrawFrame();
      if (!im.empty())
      {
        if (mbSaveResults)
        {
          std::ostringstream out;
          out << std::internal << std::setfill('0') << std::setw(5)
              << (unsigned int)(timestamp);
          cv::imwrite("2D" + out.str() + ".png", im);
        }
        cv::imshow("ORBSLAM2: Current Frame", im);
        cv::waitKey(10);
      }

      pangolin::FinishFrame();
      unique_lock<mutex> lock22(mMutexTimeStamp);

      if ((menuNext) || (menuAutoPlay))
      {
        unique_lock<mutex> locknext(mMutexNext);
        this->next = true;
        menuNext = false;
      }
      else
      {
        unique_lock<mutex> locknext(mMutexNext);
        this->next = false;
      }

      if (menuReset)
      {
        menuShowGraph = true;
        menuShowKeyFrames = false;
        menuShowPoints = true;
        menuLocalizationMode = false;
        if (bLocalizationMode)
          mpSystem->DeactivateLocalizationMode();
        bLocalizationMode = false;
        bFollow = true;
        menuFollowCamera = false;
        mpSystem->Reset();
        menuReset = false;
      }

      if (Stop())
      {
        while (isStopped())
        {
            this_thread::sleep_for(chrono::microseconds(3000));
        }
      }

      if (CheckFinish())
        break;
    }
    // MeshDrawer::Exit();
    SetFinish();
  }

  void Viewer::RequestFinish()
  {
    unique_lock<mutex> lock(mMutexFinish);
    mbFinishRequested = true;
  }

  void Viewer::RequestStop()
  {
    unique_lock<mutex> lock(mMutexStop);
    if (!mbStopped)
      mbStopRequested = true;
    this->mpMapDrawer->reset();
  }

  bool Viewer::isFinished()
  {
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinished;
  }

  bool Viewer::isStopped()
  {
    unique_lock<mutex> lock(mMutexStop);
    return mbStopped;
  }

  void Viewer::Release()
  {
    unique_lock<mutex> lock(mMutexStop);
    mbStopped = false;
  }

  void Viewer::Updatetimestamp(double t)
  {
    unique_lock<mutex> lock2(mMutexTimeStamp);
    timestamp = t;
  }

  bool Viewer::go()
  {
    unique_lock<mutex> locknext(mMutexNext);
    return this->next;
  }

  bool Viewer::Stop()
  {
    unique_lock<mutex> lock(mMutexStop);
    unique_lock<mutex> lock2(mMutexFinish);

    if (mbFinishRequested)
      return false;
    else if (mbStopRequested)
    {
      mbStopped = true;
      mbStopRequested = false;
      return true;
    }

    return false;
  }

  bool Viewer::CheckFinish()
  {
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinishRequested;
  }

  void Viewer::SetFinish()
  {
    unique_lock<mutex> lock(mMutexFinish);
    mbFinished = true;
  }
} // namespace ORB_SLAM2
