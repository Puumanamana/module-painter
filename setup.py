from setuptools import setup, find_packages

setup(name='module-painter',
      version='0.1.0',
      description='Tool for bacteriophage painting',
      long_description=open('README.rst').read(),
      url='https://github.com/Puumanamana/module-painter',
      packages=find_packages(),
      license='Apache License 2.0',
      zip_safe=False,
      entry_points={'console_scripts': ['module-painter=module_painter.main:main']},
      python_requires='>=3.6',
      test_requires=['pytest', 'pytest-cov'],
      install_requires=[
          'numpy',
          'pandas',
          'Biopython',
          'argparse',
          'bokeh>=1.4',
          'holoviews',
          'scikit-learn',
          'igraph',
          'mappy',
          'psutil',
          'seaborn'
      ],
      extras_require = {
          'full': ['pycairo']
      })
