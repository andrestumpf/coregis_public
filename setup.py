from distutils.core import setup

setup(
    name='coregis',
    version='0.1',
    packages=['rios-1.4.3.rios', 'rios-1.4.3.rios.parallel', 'rios-1.4.3.rios.riostests', 'rios-1.4.3.build.lib.rios',
              'rios-1.4.3.build.lib.rios.parallel', 'rios-1.4.3.build.lib.rios.riostests',
              'python-fmask-0.4.3.build.lib.linux-x86_64-3.5.fmask', 'python-fmask-0.4.3.fmask'],
    install_requires=[
       'gdal>=2.1.3',
       'numpy>=1.12',
       'os>=0.1>=0.1',
       'pprint>=0.1',
       'cv2',
       'matplotlib>=2.0',
       'PIL>=1.1.6',
       'shapely>=1.6'],
    url='',
    license='',
    author='André Stumpf',
    author_email='andre.stumpf@unistra.fr',
    description='A Python interface for sub-pixel co-registration of Sentinel-2 and Landsat-8 imagery'

)
