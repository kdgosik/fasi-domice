FROM gitpod/workspace-full-vnc:2022-11-15-17-00-18

# Qt5 graphics libraries for napari
RUN sudo apt-get update && \
    sudo apt-get install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools s3fs lftp && \
    sudo rm -rf /var/lib/apt/lists/*

# Napari image viewer and libraries for image processing
RUN python -m pip install --no-cache-dir \
  "napari[all]" \
  imageio-ffmpeg \
  scikit-learn \
  matplotlib \
  scikit-image \
  opencv-python \
  zarr \
  ome_types \
  basicpy \
  jupyterlab

## Nextflow
RUN /bin/bash -c 'source /home/gitpod/.sdkman/bin/sdkman-init.sh && \
  curl -s https://get.nextflow.io | bash && \
  sudo mv nextflow /usr/local/bin'
