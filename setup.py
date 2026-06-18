from glob import glob
import os
import sys
import subprocess
from setuptools import setup
import tempfile
from pybind11.setup_helpers import Pybind11Extension, build_ext
try:
    from distutils.errors import CCompilerError
except ImportError:
    from setuptools._distutils.errors import CCompilerError


__version__ = "4.0.2"

def __read__(file_name):
    return open(os.path.join(os.path.dirname(__file__), file_name)).read()


#Try to find Homebrew's libomp on macOS
def get_homebrew_libomp_prefix():
  try:
    result = subprocess.run(
      ["brew", "--prefix", "libomp"],
      capture_output=True, text=True, timeout=10)
    if result.returncode == 0:
      prefix = result.stdout.strip()
      if os.path.isdir(prefix):
        return prefix
  except (FileNotFoundError, subprocess.TimeoutExpired):
    pass
  return None


#Custom build class that provides additional checks for OpenMP
class custom_build_ext(build_ext):

  #Compile a test program to determine if C++ compiler supports OpenMP
  def check_openmp_support(self):
    os.makedirs(self.build_temp, exist_ok=True)

    #Determine platform-specific OpenMP flags
    if sys.platform == "darwin":
      libomp_prefix = get_homebrew_libomp_prefix()
      if libomp_prefix:
        print("Found Homebrew libomp at: " + libomp_prefix)
        self._omp_compile_flags = [
          "-Xclang", "-fopenmp",
          "-I" + libomp_prefix + "/include"]
        self._omp_link_flags = [
          "-L" + libomp_prefix + "/lib",
          "-Wl,-rpath," + libomp_prefix + "/lib",
          "-lomp"]
      else:
        print("Homebrew libomp not found, OpenMP will be disabled")
        return False
    else:
      self._omp_compile_flags = ["-fopenmp"]
      self._omp_link_flags = ["-fopenmp"]

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
                                        extra_postargs=self._omp_compile_flags,
                                        output_dir="/")
      except CCompilerError:
        print("Compiler does not support OpenMP")
        return False
      else:
        print("Enabling OpenMP support")

        for o in objects:
          os.remove(o)

      return True

  #Add OpenMP compiler and linker flags if necessary
  def build_extensions(self):
    use_openmp = self.check_openmp_support()

    if use_openmp:
      for ext in self.extensions:
        if not ext.extra_compile_args:
          ext.extra_compile_args = []
        ext.extra_compile_args.extend(self._omp_compile_flags)

        if not ext.extra_link_args:
          ext.extra_link_args = []
        ext.extra_link_args.extend(self._omp_link_flags)

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
    cxx_std = 17,
    language ='c++',
  ),
]


setup(
  name             = "pyfastchem",
  python_requires  = ">=3.10",
  description = "FastChem, an ultra-fast equilibrium chemistry",
  long_description=__read__('README.md'),
  long_description_content_type="text/markdown",
  author      = "Daniel Kitzmann, Joachim Stock, Brett Morris",
  license     = "GPL 3.0",
  url         = "https://github.com/NewStrangeWorlds/FastChem",
  version     = __version__,
  ext_modules = ext_modules,
  cmdclass    = {"build_ext": custom_build_ext}
)
