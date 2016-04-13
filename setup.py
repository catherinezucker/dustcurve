from setuptools import setup

setup(name='dustcurve',
      version='0.1',
      description='Refine Galactic Rotation Curve using Dust',
      url='http://github.com/phys201-sp2016/dustcurve',
      author='Catherine Zucker',
      license='GNU Public License',
      author_email='catherine.zucker@cfa.harvard.edu',
      packages=['dustcurve'],
      install_requires=[
          'numpy',
          'matplotlib',
          'emcee',
          'scipy'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
