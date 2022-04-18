from setuptools import setup, Extension, find_packages

setup(name='myspkmeans',
      version='1.0',
      install_requires=["pandas",
                        "numpy"],
      packages=find_packages(),
      description='spkmeans algorithm implementation in C',
      ext_modules=[Extension('myspkmeans', sources=['spkmeansmodule.c', 'spkmeans.c'])])
