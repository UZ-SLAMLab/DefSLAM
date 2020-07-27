echo "Uncompress vocabulary ..."

cd Vocabulary
tar -xf ORBvoc.txt.tar.gz
cd ..

echo "Configuring and building DefSLAM ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 2
