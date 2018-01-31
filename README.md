coregis: The Python implementation of <cite>[1]</cite> to improve the co-registration of 
Sentinel-2 and Landsat-8 images with a particular focus on Earth 
surface motion measurements.

[1]: Stumpf, A., Michéa, D., Malet, J.-P. (2018): Improved Co-Registration of 
Sentinel-2 and Landsat-8 Imagery for Earth Surface Motion Measurements.
Remote Sensing. 2018, 10(2), 160; doi:10.3390/rs10020160 

## Install dependencies (Ubuntu 16.04)

* Install some system libraries
```bash
sudo apt-get update
sudo apt-get install git \
                     cmake \
                     make \
                     imagemagick \
                     libimage-exiftool-perl \
                     exiv2 \
                     libgeo-proj4-perl \
                     libx11-dev \
                     cmake-curses-gui \
                     freeglut3-dev \
                     libpcre3 \
                     libpcre3-dev \
                     libtiff5 \
                     libtiff5-dev \
                     libgdal-dev
```

* Install the MicMac library
```bash
cd ~
git clone https://github.com/micmacIGN/micmac.git
cd micmac
mkdir build
cd build
cmake -DWITH_QT5=OFF -DBUILD_POISSON=OFF ../
make install -j4
echo 'export PATH="~/micmac/bin:$PATH"' >> ~/.bashrc
source .bashrc
```

* Install the Orfeo toolbox
```bash
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install otb-bin
```

* Install Anaconda
```bash
cd ~
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
bash Anaconda3-5.0.1-Linux-x86_64.sh -b -p ~/anaconda
rm Anaconda3-5.0.1-Linux-x86_64.sh
echo 'export PATH="~/anaconda/bin:$PATH"' >> ~/.bashrc 
source .bashrc
conda config --add channels conda-forge
conda update conda
```

* Create conda virtual environment with all necessary packages
```bash
conda create -n conda_coregis python=3 python-fmask opencv scipy shapely paramiko pillow paramiko matplotlib statsmodels
source activate conda_coregis
source deactivate conda_coregis
```

* Clone this repository
```bash
cd ~
git clone https://github.com/andrestumpf/coregis_public.git
```

* Compile MPIC binaries
```bash
cd ~/coregis_public/mpic
make
```

* Open it in your favourite Python IDE e.g. Pycharm
```bash
sudo apt-get install snap
sudo snap install pycharm-community --classic
```

* Make sure you use the previously setup conda environment.
In PyCharm this can be done by going to
Settings > Project > Project Interpreter > Button on the upper right >
Add local > Existing Environment > Select:
```~/anaconda/envs/conda_coregis/bin/python```




## Documentation
A more complete documentation is coming soon. For a simple use case have a look at Example.py for the moment.

# Examples #

* see Example.py
* Note that the test data is not included in the repository. 
It can be downloaded from <https://scihub.copernicus.eu>

![before coregistration](figs/regression_before_coregis.png?raw=true)
*Fig.1: Regression among master and slave image before correction.*
![after coregistration](figs/regression_after_coregis.png?raw=true)
*Fig.2: Regression among master and slave image after correction.*
# Licence #

* GPL (v3+)