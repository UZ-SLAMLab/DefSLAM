rem @echo off

SETLOCAL
set exe_dir=%~dp0..\Apps\Release\DefSLAMGT.exe
set voc_dir=%~dp0..\Vocabulary\ORBvoc.txt
set yaml_dir=%~dp0..\..\MandalaDataset\stereo0.yaml
set img_dir=\images
set time_dir=\timestamps\timestamps.txt

for %%a in (%*) do (
	echo "%%~dpna%img_dir%"
	powershell -Command %exe_dir% %voc_dir% %yaml_dir% %%~dpna%img_dir% %%~dpna%img_dir% %%~dpna%time_dir%
 )
ENDLOCAL
pause