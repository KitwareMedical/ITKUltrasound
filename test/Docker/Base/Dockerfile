FROM ipython/notebook:latest
MAINTAINER Matt McCormick <matt.mccormick@kitware.com>

RUN apt-get update && apt-get install -y \
  bash \
  build-essential \
  bzip2 \
  curl \
  git \
  libcurl4-openssl-dev \
  libssl-dev \
  make \
  ninja-build \
  python \
  python-dev \
  python-pip \
  python-numpy \
  python-scipy \
  python-nose \
  python-setuptools \
  sed \
  tar \
  vim \
  wget && \
  apt-get build-dep -y python-matplotlib

RUN pip install matplotlib

# Build and install CMake from source.
WORKDIR /usr/src
RUN git clone git://cmake.org/cmake.git CMake && \
  cd CMake && \
  git checkout v3.3.0
RUN mkdir CMake-build
WORKDIR /usr/src/CMake-build
RUN /usr/src/CMake/bootstrap \
    --parallel=$(nproc) \
    --prefix=/usr && \
  make -j$(nproc) && \
  ./bin/cmake -DCMAKE_USE_SYSTEM_CURL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_USE_OPENSSL:BOOL=ON . && \
  make install && \
  rm -rf *

WORKDIR /usr/src
RUN git clone git://itk.org/ITK.git
WORKDIR /usr/src/ITK
RUN git fetch http://review.source.kitware.com/ITK refs/changes/56/19856/5 && \
  git checkout FETCH_HEAD -b ExternalModule
RUN mkdir /usr/src/ITK-build
WORKDIR /usr/src/ITK-build
RUN cmake -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILD_TESTING:BOOL=OFF \
  -DBUILD_EXAMPLES:BOOL=OFF \
  -DBUILD_SHARED_LIBS:BOOL=ON \
  -DITK_WRAP_PYTHON:BOOL=ON \
  /usr/src/ITK && \
  ninja

RUN mkdir /usr/bin/ITKUltrasound-build

# IPython
ENV PASSWORD wavy
WORKDIR /notebooks
CMD ["/notebook.sh"]
