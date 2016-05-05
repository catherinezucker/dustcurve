from setuptools import setup

setup(name='dustcurve',
      version='0.1',
      description='Refine Galactic Rotation Curve using Dust',
      url='http://github.com/phys201-sp2016/dustcurve',
      author='Catherine Zucker',
      author_email='catherine.zucker@cfa.harvard.edu',
      packages=['dustcurve'],
      requires=['h5py', 'matplotlib', 'emcee', 'numpy', 'scipy', 'seaborn', 'pandas'],
      install_requires=[
          'h5py>=2.5.0',
          'matplotlib',
          'emcee',
          'numpy',
          'scipy',
          'seaborn',
          'pandas'
      ],
      dependency_links=['https://github.com/h5py/h5py.git'],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
      
      
