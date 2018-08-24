'''
Created on 09.11.2009

@author: philkraf
'''

from distutils.core import setup, Extension

ext = Extension('decomp._decomp',
                sources=['decomp/DECOMP.i', 'decomp/SOM.cpp', 'decomp/SOMcomponent.cpp'],
                extra_compile_args=['/openmp','/EHsc'],
                swig_opts=['-c++', '-Wextra', '-w512', '-w511', '-O', '-keyword', '-castmode'],
                )
setup(name='decomp',
      author='Philipp Kraft',
      ext_modules=[ext],
      packages=['decomp'])
