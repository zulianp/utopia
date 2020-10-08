lynx -source rawgit.com/transcode-open/apt-cyg/master/apt-cyg > apt-cyg
if ! type "$apt-cyg" > /dev/null; then
  install apt-cyg /bin
fi
apt-cyg install gcc-fortran
apt-cyg install gcc-g++
apt-cyg install gdb
apt-cyg install make
apt-cyg install cmake
apt-cyg install openmpi
apt-cyg install libopenmpi-devel
apt-cyg install libopenmpicxx1
apt-cyg install libopenblas
apt-cyg install liblapack-devel
apt-cyg install python
apt-cyg install doxygen
apt-cyg install zlib