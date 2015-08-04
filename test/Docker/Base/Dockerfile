FROM ipython/notebook:latest
MAINTAINER Matt McCormick <matt.mccormick@kitware.com>

RUN apt-get update && apt-get install -y \
  bash \
  build-essential \
  bzip2 \
  curl \
  git \
  libcurl4-openssl-dev \
  libncurses5-dev \
  libfftw3-dev \
  libssl-dev \
  make \
  ninja-build \
  python \
  python-dev \
  python-pip \
  python-matplotlib \
  python-numpy \
  python-scipy \
  python-nose \
  python-setuptools \
  sed \
  tar \
  vim \
  wget \
  libexpat1-dev \
  libhdf5-dev \
  libgdcm2-dev \
  libjpeg-dev \
  libpng12-dev \
  libtiff5-dev \
  zlib1g-dev

# Build and install CMake from source.
WORKDIR /usr/src
RUN git clone git://cmake.org/cmake.git CMake && \
  cd CMake && \
  git checkout v3.3.0 && \
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

# Build and install ITK from source.
RUN git clone git://itk.org/ITK.git && \
  cd /usr/src/ITK && \
  git fetch http://review.source.kitware.com/ITK refs/changes/56/19856/7 && \
    git checkout FETCH_HEAD -b ExternalModule && \
  mkdir /usr/src/ITK-build && \
  cd /usr/src/ITK-build && \
  cmake -G Ninja \
    -DCMAKE_BUILD_TYPE:STRING=MinSizeRel \
    -DCMAKE_INSTALL_PREFIX:PATH=/usr \
    -DBUILD_TESTING:BOOL=OFF \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
    -DITKGroup_Filtering:BOOL=ON \
    -DITKGroup_Numerics:BOOL=ON \
    -DITKGroup_IO:BOOL=ON \
    -DITK_WRAP_PYTHON:BOOL=ON \
    -DITK_USE_FFTWD:BOOL=ON \
    -DITK_USE_SYSTEM_FFTW:BOOL=ON \
    -DITK_USE_SYSTEM_SWIG:BOOL=OFF \
    -DITK_USE_SYSTEM_GDCM:BOOL=ON \
    -DITK_USE_SYSTEM_LIBRARIES:BOOL=ON \
      /usr/src/ITK && \
  ninja install && \
  cd .. && \
  rm -rf ITK ITK-build

RUN mkdir /usr/bin/ITKUltrasound-build

# IPython
ENV PASSWORD wavy
WORKDIR /notebooks
CMD ["/notebook.sh"]
