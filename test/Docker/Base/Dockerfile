FROM jupyter/scipy-notebook:latest
MAINTAINER Matt McCormick <matt.mccormick@kitware.com>

USER root
RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  bash \
  build-essential \
  bzip2 \
  curl \
  file \
  git \
  libcurl4-openssl-dev \
  libncurses5-dev \
  libfftw3-dev \
  libgl1-mesa-dev \
  libssl-dev \
  make \
  ninja-build \
  sed \
  tar \
  vim \
  wget \
  libexpat1-dev \
  libhdf5-dev \
  libjpeg-dev \
  libpng12-dev \
  libtiff5-dev \
  xserver-xorg-video-dummy \
  xserver-xorg-input-void \
  libxt-dev \
  zlib1g-dev

# Build and install CMake from source.
WORKDIR /usr/src
RUN git clone git://cmake.org/cmake.git CMake && \
  cd CMake && \
  git checkout v3.8.0 && \
  mkdir ../CMake-build && \
  cd ../CMake-build && \
  /usr/src/CMake/bootstrap \
    --parallel=$(nproc) \
    --prefix=/usr && \
  make -j$(nproc) && \
  ./bin/cmake -DCMAKE_USE_SYSTEM_CURL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_USE_OPENSSL:BOOL=ON . && \
  make install && \
  cd .. && \
  rm -rf CMake CMake-build

# Build and install KWStyle from source.
RUN git clone https://github.com/Kitware/KWStyle.git && \
  mkdir KWStyle-build && \
  cd KWStyle-build && \
  cmake \
    -G Ninja \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
      ../KWStyle && \
  ninja install && \
  cd .. && \
  rm -rf KWStyle KWStyle-build

# Build and install VTK from source
RUN git clone https://github.com/Kitware/VTK.git VTK && \
  cd VTK && git checkout v7.1.0 && cd .. && \
  mkdir VTK-build && \
  cd VTK-build && \
  cmake \
    -G Ninja \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DVTK_Group_Rendering:BOOL=OFF \
    -DVTK_WRAP_PYTHON:STRING=ON \
    -DPYTHON_EXECUTABLE:FILEPATH=/opt/conda/bin/python \
    -DVTK_PYTHON_VERSION:STRING=3 \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
      ../VTK && \
  ninja && \
  ninja install && \
  cd .. && rm -rf VTK VTK-build

# Build and install ITK from source.
# 2017-04-19 master
RUN git clone https://github.com/InsightSoftwareConsortium/ITK.git && \
  cd /usr/src/ITK && \
  git checkout 46350ec90671385d688e71aaa7058030bca6e23a && \
  mkdir /usr/src/ITK-build && \
  cd /usr/src/ITK-build && \
  cmake -G Ninja \
    -DCMAKE_BUILD_TYPE:STRING=MinSizeRel \
    -DPYTHON_EXECUTABLE:FILEPATH=/opt/conda/bin/python \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
    -DITKGroup_Filtering:BOOL=ON \
    -DITKGroup_Numerics:BOOL=ON \
    -DITKGroup_IO:BOOL=ON \
    -DModule_ITKHDF5:BOOL=ON \
    -DITK_WRAP_PYTHON:BOOL=ON \
    -DITK_USE_FFTWD:BOOL=ON \
    -DITK_USE_FFTWF:BOOL=ON \
    -DITK_USE_SYSTEM_FFTW:BOOL=ON \
    -DITK_USE_SYSTEM_SWIG:BOOL=OFF \
      /usr/src/ITK && \
  ninja && \
  cp /usr/src/ITK-build/Wrapping/Generators/Python/WrapITK.pth /opt/conda/lib/python3.5/site-packages && \
  find ./Wrapping -name '*.cpp' -delete -o -name '*.xml' -delete && \
  find . -name '*.o' -delete

RUN mkdir /usr/src/ITKUltrasound-build

USER jovyan
WORKDIR /home/jovyan
