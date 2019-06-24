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
    libassimp-dev \
    libglfw3-dev \
    libglew-dev \
    && apt-get clean \
    && apt-get autoclean \
    && apt-get autoremove


WORKDIR /first_raytracer

ARG RENDERER_VERSION=0
RUN git clone https://github.com/jammm/first_raytracer.git
RUN cd first_raytracer/first_ray \
    && mkdir deps \
    && wget https://github.com/cpp-taskflow/cpp-taskflow/archive/v2.1.0.tar.gz \
    && tar -xvf v2.1.0.tar.gz -C deps/

WORKDIR /first_raytracer/first_raytracer/first_ray

RUN g++ -O3 -std=c++17 *.cpp -o main -I./deps/cpp-taskflow-2.1.0 -lGLEW -lassimp -lglfw -lGL -pthread


### configure headless X server
COPY xorg.conf /etc/X11/xorg.conf
ENV DISPLAY :0
