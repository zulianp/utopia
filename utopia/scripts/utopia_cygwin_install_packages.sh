lynx -source rawgit.com/transcode-open/apt-cyg/master/apt-cyg > apt-cyg
if ! type "$apt-cyg" > /dev/null; then
  install apt-cyg /bin
fi
apt-cyg install gcc-fortran