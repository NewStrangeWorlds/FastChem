from glob import glob
import os
import sys
from setuptools import setup
import tempfile
from pybind11.setup_helpers import Pybind11Extension, build_ext
from distutils.dir_util import mkpath



__version__ = "2.0"


class custom_build_ext(build_ext):
  
  #Compile a test program to determine if C++ compiler supports OpenMP.
  def _check_openmp_support(self):
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
                      for (int i=0; i<10; ++i);
                    
                      return omp_get_num_threads();
                    }
                    """)
     
      srcfile.flush()
            
      try:
        objects = self.compiler.compile([srcfile.name],
                                        extra_postargs=["-fopenmp"],
                                        output_dir="/")
      except _CCompilerError:
        print("Compiler does not support OpenMP")
        use_openmp = False
      else:
        print("Enabling OpenMP support")
        use_openmp = True
          
        for o in objects:
          os.remove(o)
    
      return use_openmp

  #add OpenMP compiler and linker flags if necessary
  def build_extensions(self):
    use_openmp = self._check_openmp_support()
        
    if use_openmp:
      for ext in self.extensions:
        if not ext.extra_compile_args:
          ext.extra_compile_args = []
        ext.extra_compile_args.append('-fopenmp')
        if not ext.extra_link_args:
          ext.extra_link_args = []
        
        ext.extra_link_args.append('-fopenmp')
      
      build_ext.build_extensions(self)



ext_modules = [
  Pybind11Extension(
    "pyfastchem",
    sorted(glob("fastchem_src/*.cpp") +
           glob("python/*.cpp")),
    extra_compile_args = ['-pedantic', '-MMD'],
    extra_link_args    = ['-pedantic', '-MMD'],
    cxx_std = 11,
    language='c++',
  ),
]


setup(
  name        = "pyfastchem",
  author      = "Daniel Kitzmann, Joachim Stock, Brett Morris",
  url         = "https://github.com/exoclime/FastChem",
  version     = __version__,
  ext_modules = ext_modules,
  cmdclass    = {"build_ext": custom_build_ext}
)
