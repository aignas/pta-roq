from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy import get_include

ext_modules = [
    Extension("signal_model",  ["signal_model.pyx"],  include_dirs=[get_include()]),
    Extension("reduced_basis", ["reduced_basis.pyx"], include_dirs=[get_include()], extra_compile_args=['-fopenmp'], extra_link_args=['-fopenmp'])
    ]

setup(
  name = 'PTA signal modelling',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
