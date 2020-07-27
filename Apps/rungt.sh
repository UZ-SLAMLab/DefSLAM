# Run camera
#./DefSLAM /home/jose/DeformableSLAM/Vocabulary/ORBvoc.txt /home/jose/DeformableSLAM/calibration_files/logitechc922.yaml
# Run one video without ground truth
#./DefSLAM /home/jose/DeformableSLAM/Vocabulary/ORBvoc.txt /media/jose/NuevoVol/videosDataset/sequence_heart/hamlyn.yaml /media/jose/NuevoVol/videosDataset/f5phantom/f5_dynamic_deint_L.avi

# Groundtruth depth image
#./DefSLAMGTCT /home/jose/DeformableSLAM/Vocabulary/ORBvoc.txt /media/jose/NuevoVol/videosDataset/f5phantom/hamlyn.yaml /media/jose/NuevoVol/videosDataset/f5phantom/f5_dynamic_deint_L.avi /media/jose/NuevoVol/videosDataset/f5phantom/f5/heartDepthMap_

# Groundtruth stereo
./DefSLAMGT /home/jose/DeformableSLAM/Vocabulary/ORBvoc.txt /media/jose/NuevoVol/videosDataset/Jose/stereo3.yaml /media/jose/NuevoVol/videosDataset/Jose/Mandala3/images /media/jose/NuevoVol/videosDataset/Jose/Mandala3/images /media/jose/NuevoVol/videosDataset/Jose/Mandala3/timestamps/timestamps.txt
