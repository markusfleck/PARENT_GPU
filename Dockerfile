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

FROM nvidia/cuda:10.1-devel

RUN apt-get update

# add tools here which you want to use
RUN apt-get install -y nano vim git

RUN mkdir /PARENT_GPU
COPY ./ /PARENT_GPU

RUN apt-get -y install wget cmake g++ sudo

RUN mkdir /gromacs-2018.8_build/
RUN cd /gromacs-2018.8_build/; wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-2018.8.tar.gz; tar xvfz gromacs-2018.8.tar.gz; mkdir -p gromacs-2018.8/build
RUN cd /gromacs-2018.8_build/gromacs-2018.8/build; export CC=`which gcc-8`; cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DCMAKE_INSTALL_PREFIX=/gromacs-2018.8
RUN cd /gromacs-2018.8_build/gromacs-2018.8/build; export CC=`which gcc-8`; make -j 4; make install

RUN useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo
RUN chown -R docker:docker /PARENT_GPU;
USER docker

ENV CPLUS_INCLUDE_PATH "/gromacs-2018.8/include:$CPLUS_INCLUDE_PATH"
ENV LIBRARY_PATH "/gromacs-2018.8/lib:$LIBRARY_PATH"
ENV LD_LIBRARY_PATH "/gromacs-2018.8/lib:$LD_LIBRARY_PATH"

RUN cd /PARENT_GPU; make clean; make; 

CMD cd /PARENT_GPU; echo Run \"make checks\" to test the installation. Run \"make clean\" in case of error.;/bin/bash



