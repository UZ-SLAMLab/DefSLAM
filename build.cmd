IF NOT DEFINED VCPKG_ROOT SET VCPKG_ROOT="C:\git\vcpkg"
IF NOT DEFINED VCPKG_DEFAULT_TRIPLET SET VCPKG_DEFAULT_TRIPLET=x64-windows
@del /S/Q build_win
@rmdir /S/Q build_win
@mkdir build_win
@cd build_win
@call cmake -G"Visual Studio 16 2019" -DCMAKE_TOOLCHAIN_FILE="%VCPKG_ROOT%/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=%VCPKG_DEFAULT_TRIPLET% ..
@cmake --build . --config Release 
@cd ..