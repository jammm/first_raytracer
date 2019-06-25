FROM ubuntu:18.04

### headless X server dependencies
RUN apt-get update \
    && apt-get install -y \
    libx11-dev \
    libxxf86vm-dev \
    x11-xserver-utils \
    x11proto-xf86vidmode-dev \
    x11vnc \
    xpra \
    xserver-xorg-video-dummy \
    && apt-get clean \
    && apt-get autoclean \
    && apt-get autoremove

RUN apt-get update \
    && apt-get install -y \
    git \
    build-essential \
    wget \
    xorg-dev \
    cmake \
    libassimp-dev \
    libglew-dev \
    && apt-get clean \
    && apt-get autoclean \
    && apt-get autoremove

WORKDIR /first_raytracer

RUN git clone https://github.com/jammm/first_raytracer.git

RUN cd first_raytracer/first_ray \
    && wget https://github.com/glfw/glfw/archive/3.3.tar.gz \
    && mkdir deps \
    && tar -xvf 3.3.tar.gz -C deps/ \
    && cd deps/glfw-3.3 \
    && cmake -DBUILD_SHARED_LIBS=ON . \
    && make \ 
    && make install \
    && cd ../ \
    && rm -rf ../3.3.tar.gz glfw-3.3

RUN cd first_raytracer/first_ray \
    && wget https://github.com/cpp-taskflow/cpp-taskflow/archive/v2.1.0.tar.gz \
    && tar -xvf v2.1.0.tar.gz -C deps/

WORKDIR /first_raytracer/first_raytracer/first_ray

RUN g++ -O3 -std=c++17 *.cpp -o main -I./deps/cpp-taskflow-2.1.0 -lGLEW -lassimp -lglfw -lGL -pthread


### configure headless X server
COPY xorg.conf /etc/X11/xorg.conf
ENV DISPLAY :0
COPY entrypoint.sh /

ENTRYPOINT "/entrypoint.sh"
