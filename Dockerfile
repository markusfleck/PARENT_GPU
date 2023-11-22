# You need a CUDA aware docker installation, see 
# https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html
# Build image via "docker build -t parent_gpu .".
# Run the image via "docker run -it -v $(pwd):/PARENT_GPU --name parent_gpu --gpus all parent_gpu".
# The -v flag maps your current host directory directory to the /PARENT_GPU 
# directory inside your container. Make sure you don't have binaries (obj and exec folders)
# compiled on the host system or another container in the host folder ("make clean" helps).
# After closing with "exit", reaccess your container
# via "docker container restart parent_gpu; docker container attach parent_gpu".
# sudo password for user "docker" inside container is just "docker".

FROM nvidia/cuda:12.3.0-devel-ubuntu22.04

RUN apt-get update

RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata

# add tools here which you want to use
RUN apt-get install -y nano vim git

RUN mkdir /PARENT_GPU
COPY ./ /PARENT_GPU

RUN apt-get -y install wget cmake g++ sudo
RUN apt-get -y install libgromacs-dev

RUN useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo
RUN chown -R docker:docker /PARENT_GPU;
USER docker

RUN cd /PARENT_GPU; make; 
CMD cd /PARENT_GPU; echo Run \"make checks\" to test the installation. Run \"make clean\" in case of error.;/bin/bash

