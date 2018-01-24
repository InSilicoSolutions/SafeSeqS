from setuptools import setup;
from setuptools import find_packages;

def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='safeseqs',
      version='0.91',
      description='Safe-SequencingSystem (SafeSeqs) data processing pipeline',
      long_description=readme(),
      url='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3111315/',
      author='Kenneth Kinzler',
      author_email='kinzlke@jhmi.edu',
      license='',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          'scipy'
      ],
      classifiers=[
      'Development Status :: 5 - Production/Stable',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
