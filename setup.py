from setuptools import setup, find_packages

setup(name='modular_painter',
      version='0.0',
      description='Tool for bacteriophage painting',
      long_description=open('README.rst').read(),
      url='https://github.com/Puumanamana/modular_painting',
      packages=find_packages(),
      license='Apache License 2.0',
      zip_safe=False,
      python_requires='>=3.6',
      test_requires=['pytest', 'pytest-cov'],
      install_requires=[
          'numpy',
          'pandas',
          'Biopython',
          'argparse',
          'bokeh>=1.4',
          'mappy'
      ])
