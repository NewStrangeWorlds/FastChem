from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


__version__ = "2.0"

ext_modules = [
    Pybind11Extension(
        "pyfastchem",
        sorted(glob("model_src/*.cpp") +
               glob("fastchem_src/*.cpp") +
               glob("python/*.cpp")),
    ),
]

setup(
    name="pyfastchem",
    author="Daniel Kitzmann, Joachim Stock, Brett Morris",
    url="https://github.com/exoclime/FastChem",
    version=__version__,
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext}
)
