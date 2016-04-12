from setuptools import setup

setup(name='dustcurve',
      version='0.1',
      description='Create Dust-Refined Galactic Rotation Curve',
      url='http://github.com/p201-sp2016/Refining_Kinematic_Distance_Estimates',
      author='Catherine Zucker',
      author_email='catherine.zucker@cfa.harvard.edu,
      license='None',
      packages=['dustcurve'],
      install_requires=[
          'numpy',
      ],
      zip_safe=False)
