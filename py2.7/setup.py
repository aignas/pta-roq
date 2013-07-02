from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension("signal_model",  ["signal_model.pyx"]),
    Extension("reduced_basis", ["reduced_basis.pyx"], extra_compile_args=['-fopenmp'], extra_link_args=['-fopenmp'])
    ]

setup(
  name = 'PTA signal modelling',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
