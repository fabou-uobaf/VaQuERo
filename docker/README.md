# Running VaQuERo in a container

VaQuERo requires a  number of R packages properly installed, which again make use of system libraries. To facilitate installation of the former a package installation script is provided which can be invoked by

``Rscript --vanilla VaQuERo/R_package_dependency_install.r``

To provide a the proper system libraries and support reproducibility, we provide a Docker image equipped with Ubuntu 20.04 image augmented with all system dependencies.

## Run analysis within Docker container

For convenience an pre-processed image is available under https://hub.docker.com/r/wasa000/vaquero and can be utilized directly with the following command.

``sudo docker run -it -v $(pwd):/home/user/current_analysis --name VaQuERo_Container wasa000/vaquero``

To analyze your own data take care that your data are inside the current working directory where you invoke the ``docker run`` run command or adapt the command substituting the ``$(pwd)`` with the path to the appropriate data directory.

Within the container the latest version of VaQuERo can be retrieved by

``git clone https://github.com/fabou-uobaf/VaQuERo.git``

And a test run with the test data provided in the git repo can be performed by

``Rscript  --vanilla VaQuERo/scripts/VaQuERo.r --debug TRUE``

## Create own image file

A suitable docker image can be created according the specification in the attached Dockerfile by executing:

``sudo docker build -t VaQuERo_Container VaQuERo/docker/Dockerfile``
