from glob import glob
import os, zipfile
import sys
from setuptools import setup
import tempfile
from pybind11.setup_helpers import Pybind11Extension, build_ext
from distutils.dir_util import mkpath
from distutils.errors import CCompilerError, CompileError
from distutils import sysconfig
from urllib.request import urlretrieve


__version__ = "2.2"


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


  def check_boost_math(self, git_tag):
    print("Checking for Boost Math library...")
    
    dep_dir = '_deps/'
    boost_lib_dir = dep_dir + 'boost_math-src/'
    
    
    if not os.path.isdir(boost_lib_dir):
      print('Boost Math library not found. Downloading version', git_tag, 'from GitHub.')
      if not os.path.isdir(dep_dir):
        os.mkdir(dep_dir)
        
      zip_filename = 'boost_math.zip'
      git_dirname  = 'math-' + git_tag
      
      try:
        urlretrieve('https://github.com/boostorg/math/archive/' + git_tag + '.zip', zip_filename)
    
        zip_file = zipfile.ZipFile(zip_filename, 'r')
        zip_file.extractall(path=dep_dir)
        zip_file.close()
    
        os.rename(dep_dir + git_dirname, boost_lib_dir)
        os.remove(zip_filename)
      except Exception as e:
        raise CompileError(str(e) + "\nError downloading Boost math library, aborting...\n")
    
    return None


  #Add OpenMP compiler and linker flags if necessary
  def build_extensions(self):
    use_openmp = self.check_openmp_support()
    
    boost_math_git_tag = 'ed01dae24893bb69c02c6d599acc74bdb8f46bda'
    self.check_boost_math(boost_math_git_tag)

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
           glob("python/fastchem_python_wrapper.cpp")),
    define_macros = [('_SETUP_PY', '1')],
    include_dirs = ['_deps/boost_math-src/include/'],
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
