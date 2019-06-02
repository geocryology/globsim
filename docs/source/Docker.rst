.. _docker:

Running GlobSim with Docker
============================

As an alternative to building GlobSim from scratch, you can run it in a Docker container and save yourself hours of compiling ESMF.

Install Docker
--------------
Install Docker from the `Docker website <https://www.docker.com/get-started>`_


Open GlobSim Docker container
-----------------------------
Next, open a docker command-line window and download the globsim docker image::

    docker pull geocryology/globsim:v1_0_5_GMDD
    
Then, start the docker container! You'll want to attach it to a folder on your computer using the `-v` flag, replacing <directory> with a directory on your computer::

    docker run --rm -v <directory>:/data: -t -i GlobSim
    
.. note:: Windows users: Docker is only able to attach directories that are listed in VirtualBox as shared directories. This defaults to 'C:/Users' (entered as '/c/users' in the docker command line). For any directories outside of this path, you must `manually add shared folders <http://support.divio.com/local-development/docker/how-to-use-a-directory-outside-cusers-with-docker-toolboxdocker-for-windows>`_ through the VirtualBox GUI.  
.. warning:: Docker containers do not store data on them once you quit!  Make sure you keep the docker container open until the downloading and processing are completed, and move any important data onto the host machine using the attached volume!

Download GlobSim (optional)
---------------------------
The Docker container comes pre-loaded with globsim, but if you want the latest version you can download it:

From github::
  
    cd /opt
    git clone https://geocryology/globsim
    git checkout working
    
With pip (not yet)::
 
    pip3 install globsim
    
From a tar.gz file in your attached volume::

    cp /data/globsim.tar.gz /opt/
    cd /opt
    tar -xvf globsim.tar.gz

.. note:: For the examples to work, you must have the globsim directory in the docker /opt directory!

Copy credentials files
----------------------
Obtain the necessary :ref:`credentials`. files and copy them into the /root folder of the docker image using the attached volume (if you were running this elsewhere you would copy them into your home folder)

Run GlobSim using the example directories
-----------------------------------------
You're now ready to run GlobSim! Try the :ref:`example`

