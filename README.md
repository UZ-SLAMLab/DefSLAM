# 1. DefSLAM
**Authors:** [Jose Lamarca](http://webdiis.unizar.es/~jlamarca/), Shaifali Parashar, Adrien Bartoli and [J. M. M. Montiel](http://webdiis.unizar.es/~josemari/) 

DefSLAM is a real-time deformable SLAM library for **Monocular** cameras that computes the camera trajectory and a sparse 3D reconstruction in a deformable environment.

[![](https://img.youtube.com/vi/6mmhD2_t6Gs/0.jpg)](https://www.youtube.com/watch?v=6mmhD2_t6Gs&t=3s)

Find more videos of experiments here: https://www.youtube.com/playlist?list=PLKBuKNhAV30SlKGJ9eaMlAExdWRypUy-K
### Related Publications:
[Lamarca, Jose, et al. "DefSLAM: Tracking and Mapping of Deforming Scenes from Monocular Sequences." arXiv preprint arXiv:1908.08918 (2019).](https://arxiv.org/abs/1908.08918)
```
@article{lamarca2019defslam,
  title={DefSLAM: Tracking and Mapping of Deforming Scenes from Monocular Sequences},
  author={Lamarca, Jose and Parashar, Shaifali and Bartoli, Adrien and Montiel, JMM},
  journal={arXiv preprint arXiv:1908.08918},
  year={2019}
}
```

[Lamarca J., Montiel J.M.M. (2019) Camera Tracking for SLAM in Deformable Maps. In: Leal-Taixé L., Roth S. (eds) Computer Vision – ECCV 2018 Workshops. ECCV 2018. Lecture Notes in Computer Science, vol 11129. Springer, Cham](https://openaccess.thecvf.com/content_ECCVW_2018/papers/11129/Lamarca_Camera_Tracking_for_SLAM_in_Deformable_Maps_ECCVW_2018_paper.pdf)
```
@inproceedings{lamarca2018camera,
  title={Camera tracking for SLAM in deformable maps},
  author={Lamarca, Jose and Montiel, JMM},
  booktitle={Proceedings of the European Conference on Computer Vision (ECCV)},
  year={2018}
}
```
# 2. Prerequisites
We have tested the library in **16.04** and **18.04**, but it should be easy to compile in other platforms. A powerful computer (e.g. i7) will ensure real-time performance and provide more stable and accurate results.

## Pangolin
We use [Pangolin](https://github.com/stevenlovegrove/Pangolin) for visualization and user interface. Dowload and install instructions can be found at: https://github.com/stevenlovegrove/Pangolin.

## OpenCV
We use [OpenCV](http://opencv.org) to manipulate images and features. Dowload and install instructions can be found at: http://opencv.org. **Required 4.0.0**.

## Eigen3
Required by g2o (see below). Download and install instructions can be found at: http://eigen.tuxfamily.org. **Required at least 3.1.0**.

# Ceres library:
We use [Ceres](http://opencv.org) to optimize warp and to run the NRSfM.

## PCL
It is used just for ground truths. It is not critic for the program. We use [PCL](https://pointclouds.org/downloads/) mainly for groundtruths.

## DBoW2 and g2o and ORBSLAM(Included in Thirdparty folder)
We use modified versions of the [DBoW2](https://github.com/dorian3d/DBoW2) library to perform place recognition and [g2o](https://github.com/RainerKuemmerle/g2o) library to perform non-linear optimizations. Both modified libraries (which are BSD) are included in the *Thirdparty* folder. Likewise we include a modified version of ORBSLAM. BOW is only used for ORBSLAM_2.


# 3. License

DefSLAM is released under a [GPLv3 license](https://github.com/jlamraca/DefSLAM/blob/master/License-gpl.txt). For a list of all code/library dependencies (and associated licenses), please see [Dependencies.md](https://github.com/jlamarca/DefSLAM/blob/master/Dependencies.md).


# 4. Building DefSLAM library and examples

Clone the repository:
```
git clone https://github.com/JoseLamarca/DefSLAM.git
```

We provide a script `build.sh` to build the *DefSLAM* including. Please make sure you have installed all required dependencies (see section 2). Execute:
```
cd DefSLAM
chmod +x build.sh
./build.sh
```

This will create **libDefSLAM.so**  at *lib* folder and the executables **simplestereo** **simpleCT** and **simple** in *Apps* folder.

# 5. Datasets
## Mandala dataset
Download the sequences with the link: 
[Mandala Dataset](https://drive.google.com/file/d/1i3i2f3Ph22DfZ6AfXKjPRb8WrGNw_41C/view?usp=sharing)

It contains the five sequences presented in the paper. The mandala deformation becomes more challenging with the number of dataset.

The dataset is organized as follows:
```
-- yalm file
-- Mandala(i)
   -- images
   -- times
```
## Hamlyn dataset  
1. Download the dataset from the webpage [Hamlyn](http://hamlyn.doc.ic.ac.uk/vision/). 

2. To downloaded in the format run in our paper [Hamlyn not official](https://drive.google.com/file/d/1fG6RNQlBqwDKZSRUuiEuoMbVosMo5u8J/view?usp=sharing)



## Run Example
Examples of scripts in Apps/rungt.sh

1. Process a sequence. 
```
./DefSLAM <ORBfile> <yaml.file> <ORBfile> <image folder>
```
If you run it without <image folder>. It takes the camera in 0. 

2. Sequences with depth image for ground truth. (Used for CT phantom dataset)
```
./DefSLAMCTGT <ORBfile> <yaml.file> <video.avi> <folder with pattern>
```

3. Sequences with stereo for ground truth.
```
./DefSLAMGT <ORBfile> <yaml.file> <left_image folder> <right_image folder> <timestamps file>
```

Script to generate this folders with this format avaiable. Ask me through <jlamarca@unizar.es>

## Code Organization
```
-- Root folder
  -- Apps. Folder with executables.

  -- Vocabulary. Folder that contains the vocabulary for the BoW.

  -- ThirdParty
    -- BBS. Bspline library adapted for C++ with eigen.
    -- DBoW2. Bag Of word library only used for ORBSLAM_2. Loop closin and relocalization.
    -- g2o. Optimization library for deformation tracking and pose estimation.
    -- ORBSLAM_2. Base clases of the rigid SLAM library modified to work in deforming 
		  environments.

  -- Modules. DefSLAM modules, it includes all the modules needed for the library.
    -- Common. It includes the keyframes, map, map points and Systems.
    -- Mapping.
    -- GroundTruth
    -- Matching
    -- Settings. Check advance settings.
    -- Template. Template tools
    -- ToolsPCL. Tools from PCL library used mainly for the ground truth.
    -- Tracking. 
    -- Viewer.
```
## Advance settings for and ORBSLAM:
To repeat the experiments and save the results you can activate the flag Viewer.SaveResults in the yalm file. The system is running by default in parallel. To change it, you must change the file set_MAC.h and comment the line #define PARALLEL.

To run ORBSLAM, you can uncomment the line 5 in set_MAC.h.

The file CC_MAC.h contain the parameters for the cross correlation matching used in the grundtruth. It is adapted for a better performance in the Mandala dataset, but it should be adapted for each dataset.

To set the number of nodes of the mesh (Initially 10x10), change line 63 and 64 of Template/TriangularMesh.cc

To set the number of nodes for the BBSpline Thirdparty/BBS/bbs_MAC.h
