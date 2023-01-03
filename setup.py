from glob import glob
import os
import sys
from setuptools import setup
import tempfile
from pybind11.setup_helpers import Pybind11Extension, build_ext
from distutils.dir_util import mkpath
from distutils.errors import CCompilerError
from distutils import sysconfig


__version__ = "3.0"


#Custom build class that provides additional checks for OpenMP
class custom_build_ext(build_ext):
  
  #Compile a test program to determine if C++ compiler supports OpenMP
  def check_openmp_support(self):
    mkpath(self.build_temp)
        
    with tempfile.NamedTemporaryFile(mode='w',
                                     dir=self.build_temp,
                                     prefix='openmptest',
                                     suffix='.cpp') as srcfile:
            
      print("Checking if compiler supports OpenMP")
            
      srcfile.write("""
                    #include <omp.h>
                   
                    int omp_test() 
                    {
                      #pragma omp parallel for
                      for (int i=0; i<5; ++i);
                    
                      return omp_get_num_threads();
                    }
                    """)
     
      srcfile.flush()
            
      try:
        objects = self.compiler.compile([srcfile.name],
                                        extra_postargs=["-fopenmp"],
                                        output_dir="/")
      except CCompilerError:
        print("Compiler does not support OpenMP")
        use_openmp = False
      else:
        print("Enabling OpenMP support")
        use_openmp = True
          
        for o in objects:
          os.remove(o)
    
      return use_openmp

  #Add OpenMP compiler and linker flags if necessary
  def build_extensions(self):
    use_openmp = self.check_openmp_support()
        
    if use_openmp:
      for ext in self.extensions:
        if not ext.extra_compile_args:
          ext.extra_compile_args = []
        ext.extra_compile_args.append('-fopenmp')
        
        if not ext.extra_link_args:
          ext.extra_link_args = []
        ext.extra_link_args.append('-fopenmp')
    
    #Call the build function from the parent class
    build_ext.build_extensions(self)



ext_modules = [
  Pybind11Extension(
    "pyfastchem",
    sorted(glob("fastchem_src/*.cpp") +
           glob("fastchem_src/elements/*.cpp") +
           glob("fastchem_src/gas_phase/*.cpp") +
           glob("fastchem_src/condensed_phase/*.cpp") +
           glob("python/fastchem_python_wrapper.cpp")),
    define_macros = [('_SETUP_PY', '1')],
    cxx_std = 11,
    language ='c++',
  ),
]


setup(
  name        = "pyfastchem",
  description = "FastChem, an ultra-fast equilibrium chemistry",
  author      = "Daniel Kitzmann, Joachim Stock, Brett Morris",
  license     = "GPL 3.0",
  url         = "https://github.com/exoclime/FastChem",
  version     = __version__,
  ext_modules = ext_modules,
  cmdclass    = {"build_ext": custom_build_ext}
)
