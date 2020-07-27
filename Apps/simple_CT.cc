#include <System.h>
#include <Viewer.h>
#include <opencv2/core/core.hpp>
#include <unistd.h>

void loadCT(string, cv::Mat &);

int main(int argc, char **argv)
{
  if (argc != 5)
  {
    cerr << endl
         << "Usage: ./DefSLAMCTGT ORBvocabulary calibrationFile "
            "video CTfiles(heartDepthMap_ from hamlyn)"
         << endl;
    return 1;
  }

  cv::VideoCapture cap(argv[3]); // open the default camera

  if (!cap.isOpened()) // check if we succeeded
    return -1;

  string arg = argv[2]; //../../calibration_files/logitechc922.yaml";
  string arg2 = argv[1];
  string CTs = argv[4];

  cv::FileStorage fsSettings(arg, cv::FileStorage::READ);
  if (!fsSettings.isOpened())
  {
    cerr << "ERROR: Wrong path to settings" << endl;
    return -1;
  }

  cv::Mat K_l, K_r, P_l, P_r, R_l, R_r, D_l, D_r;
  fsSettings["LEFT.K"] >> K_l;
  fsSettings["RIGHT.K"] >> K_r;
  fsSettings["LEFT.P"] >> P_l;
  fsSettings["RIGHT.P"] >> P_r;
  fsSettings["LEFT.R"] >> R_l;
  fsSettings["RIGHT.R"] >> R_r;
  fsSettings["LEFT.D"] >> D_l;
  fsSettings["RIGHT.D"] >> D_r;

  int rows_l = fsSettings["LEFT.height"];
  int cols_l = fsSettings["LEFT.width"];
  int rows_r = fsSettings["RIGHT.height"];
  int cols_r = fsSettings["RIGHT.width"];

  if (K_l.empty() || K_r.empty() || P_l.empty() || P_r.empty() || R_l.empty() ||
      R_r.empty() || D_l.empty() || D_r.empty() || rows_l == 0 || rows_r == 0 ||
      cols_l == 0 || cols_r == 0)
  {
    cerr << "ERROR: Calibration parameters to rectify stereo are missing!"
         << endl;
    return -1;
  }

  cv::Mat M1l, M2l;
  cv::initUndistortRectifyMap(K_l, D_l, R_l, P_l.rowRange(0, 3).colRange(0, 3),
                              cv::Size(cols_l, rows_l), CV_32F, M1l, M2l);

  uint i(0);
  defSLAM::System SLAM(arg2, arg, true);
  while (true)
  {
    cv::Mat imLeft;
    cap >> imLeft;
    std::cout << i << " " << (double(i) / 25.0 + 0.466667) * 30 << " "
              << int((double(i) / 25.0 + 0.466667) * 30) << " "
              << int((double(i) / 25.0 + 0.466667) * 30) % 20 << std::endl;
    uint GTNum = int((double(i) / 25.0 + 0.466667) * 30) % 20;
    cv::Mat CTdepth(imLeft.rows, imLeft.cols, CV_32FC1);
    string CTGTname = CTs + to_string(GTNum) + ".txt";
    loadCT(CTGTname, CTdepth);

    if (imLeft.empty())
      break;
    if (imLeft.channels() == 1)
    {
      cv::cvtColor(imLeft, imLeft, cv::COLOR_GRAY2RGB);
    }
    cv::Mat imLeftRect;
    cv::remap(imLeft, imLeftRect, M1l, M2l, cv::INTER_LINEAR);
    cv::remap(CTdepth, CTdepth, M1l, M2l, cv::INTER_LINEAR);
    cv::Mat _mask(imLeftRect.rows, imLeftRect.cols, CV_8UC1, cv::Scalar(255));

    SLAM.TrackMonocularCTGT(imLeftRect, CTdepth, i, _mask);
    i++;
  }

  SLAM.Shutdown();

  return 0;
}

void loadCT(string CTName, cv::Mat &depth)
{
  ifstream infile(CTName);
  int rows = depth.rows;
  int cols = depth.cols;
  if (infile.is_open())
    for (int i = 0; i < rows * cols; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        double f(0.0);
        infile >> f;
        if (j == 2)
        {
          depth.at<float>(int(i / cols), int((i % cols))) = f;
        }
      }
    }
}
