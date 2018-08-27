'''
Created on 09.11.2009

@author: philkraf
'''

from __future__ import print_function, division

import io
import re

from distutils.sysconfig import customize_compiler
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class DecompBuildExt(build_ext):
    """
    Custom build class to get rid of the -Wstrict-prototypes warning
    source: https://stackoverflow.com/a/36293331/3032680

    Removes also CLASS_swigregister clutter in decomp.py
    """
    @staticmethod
    def clean_swigregister(decomp_py):
        """
        SWIG creates ugly CLASS_swigregister functions. We remove them.
        """
        rp_call = re.compile(r'^(\w*?)_swigregister\((\w*?)\)', re.MULTILINE)
        decomp_py, n = rp_call.subn('# \\1 end', decomp_py)
        print(n, 'CLASS_swigregister(CLASS) lines deleted')
        rp_def = re.compile(r'^(\w*?)_swigregister = _decomp\.(\w*?)_swigregister', re.MULTILINE)
        decomp_py, n = rp_def.subn('_decomp.\\1_swigregister(\\1)', decomp_py)
        print(n, 'CLASS_swigregister = _decomp... -> _decomp.CLASS_swigregister(CLASS)')
        return decomp_py

    @staticmethod
    def clean_static_methods(decomp_py):
        """SWIG creates still static methods as free functions (extra) to ensure Py2.2 compatibility
        We don't want that in 2018
        """
        # Find class names and free functions
        classes = re.findall(r'class\s(\w*?)\(.*?\):', decomp_py, re.MULTILINE)
        funcs = re.findall(r'^def (\w*?)\(\*args, \*\*kwargs\):$', decomp_py, flags=re.MULTILINE)

        # Find old style static methods (def CLASS_method():)
        static_methods = [f for f in funcs if [c for c in classes if f.startswith(c)]]

        count = 0
        for sm in static_methods:
            decomp_py, n = re.subn(r'^def {}.*?:.*?return.*?$'.format(sm),
                                     '\n\n', decomp_py, flags=re.MULTILINE + re.DOTALL)
            count += n
        print(count, 'old style static methods removed from', len(classes), 'classes')
        return decomp_py

    def build_extensions(self):
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        build_ext.build_extensions(self)

        decomp_py = open('decomp/decomp.py').read()
        # decomp_py = self.clean_swigregister(decomp_py)
        decomp_py = self.clean_static_methods(decomp_py)

        open('decomp/decomp.py', 'w').write(decomp_py)


ext = Extension('decomp._decomp',
                sources=['decomp/decomp.i', 'decomp/SOM.cpp', 'decomp/SOMcomponent.cpp'],
                swig_opts=['-c++', '-Wextra', '-w512', '-w511', '-O', '-keyword', '-castmode'],
                )


def get_version():
    for line in open('decomp/__init__.py'):
        if line.strip().startswith('__version__'):
            return line.split('=')[-1].strip().strip("'")


if __name__ == '__main__':

    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C++',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]

    description = 'DECOMP: a carbon pool model for soil organic matter in the forest floor'
    long_description = io.open('README.rst', encoding='utf-8').read()

    setup(name='decomp',
          version=get_version(),
          license='MIT',
          ext_modules=[ext],
          packages=['decomp'],
          python_requires='>=3.5',
          keywords='decomposition soil litter',
          author='Philipp Kraft',
          author_email="philipp.kraft@umwelt.uni-giessen.de",
          url="https://www.uni-giessen.de/hydro/download",
          description=description,
          long_description=long_description,
          classifiers=classifiers,
          cmdclass=dict(build_ext=DecompBuildExt)
          )
