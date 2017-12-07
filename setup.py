from setuptools import setup;
from setuptools import find_packages;

def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='SafeSeqS',
      version='1.1',
      description='SafeSeqs',
      long_description=readme(),
      url='http://jhu.edu',
      author='Kenneth Kinzler',
      author_email='kinzlke@jhmi.edu',
      license='',
      packages=find_packages(),
      install_requires=[
          'scipy',
          'pypiwin32'
      ],
      classifiers=[
      'Development Status :: 5 - Production/Stable',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
