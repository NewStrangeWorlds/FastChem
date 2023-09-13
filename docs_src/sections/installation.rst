Installation
============


``FastChem`` can be installed in two different ways: either using
``CMake`` or by calling a Python setup function ``setup.py``. The former
will install the ``C++`` stand-alone executable and optionally the
Python module, while the latter one will only provide the ``pyFastChem``
Python module. The Python module created by ``CMake`` will only be
available locally in the ``python`` directory, while the one produced by
``setup.py`` will be integrated in your standard Python library and,
thus, work as a normal Python package. Additionally, we also support a
Python installation via ``PyPI``, the Python Package Index.

Obtaining the code
~~~~~~~~~~~~~~~~~~

``FastChem`` is hosted on the Exoclime GitHub page:
https://github.com/exoclime/fastchem. If ``git`` is available on a
computer, the repository can be simply cloned with

.. code:: bash

   git clone https://github.com/exoclime/fastchem

Prerequisites
~~~~~~~~~~~~~

``FastChem`` is written in ``C++``. It uses features of the ``C++11``
standard and, therefore, requires a compiler that implements this
standard. We also provide an optional Python interface, allowing
``FastChem`` to be called directly from within a Python script. The
interface is based on the Python package ``PyBind11``.

Prerequisites for installation via ``CMake``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The complete list of prerequisites for a basic ``CMake`` installation
is:

-  a ``C++`` compiler (e.g. ``g++`` or ``Clang`` on MacOS)

-  ``CMake``, at least version 3.10

The ``C++`` compiler will be detected by the CMake script when it
generates the makefiles. For some of its optional components
``FastChem`` will need:

-  an ``OpenMP`` library (to run ``FastChem`` in parallel)

-  a Python 3.x interpreter (for the Python interface)

Prerequisites for Python installation via ``setup.py`` or ``PyPi`` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An installation of ``pyFastChem`` with the ``setup.py`` script
or\ ``PyPI`` requires

-  a Python 3.x interpreter

-  a ``C++`` compiler (e.g. ``g++`` or ``Clang`` on MacOS)

-  an ``OpenMP`` library (optional, required to run ``FastChem`` in
   parallel)

-  pip (when using ``PyPI``)

as well as the following Python modules:

-  PyBind11

-  setuptools

-  distutils

-  glob

-  tempfile

Supported C++ compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The compilation of ``FastChem`` has been tested on a variety of
different compilers and platforms. In particular, it was verified that
``FastChem`` can be compiled with:

-  GCC 7.5 or newer

-  Clang 12.0 (including Apple’s Clang 12.0)

Since ``FastChem`` just uses plain ``C++`` without any external library,
any compiler that supports the ``C++11`` standard should be able to
compile the code successfully.

PyBind11 Library
^^^^^^^^^^^^^^^^

For its Python interface, ``FastChem`` requires the ``PyBind11`` library
that translates the Python calls into ``C++``. While ``PyBind11`` can in
theory be installed via ``pip``, ``conda``, or ``brew`` (on MacOS),
``CMake`` isn’t always able to properly locate the installed library.

For the installation via ``CMake``, we therefore chose to include
``PyBind11`` as a submodule in the ``FastChem`` repository. ``CMake``
will take header files and Python scripts provided by the submodule to
create the ``PyFastChem`` module. No separate compilation or
installation of ``PyBind11`` is required. During the setup stage,
``CMake`` will download the ``PyBind11`` library automatically. This
code will be placed into a separate ``_deps`` folder.

If you choose to install ``pyFastChem`` via the ``setup.py`` function,
then the ``PyBind11`` library has to already present in your local
Python installation.

.. _sec:install_config:

Configuration and compilation with CMake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Before ``FastChem`` can be compiled, ``CMake`` is required to
  configure the compilation files, locate libraries, and write the
  makefiles that will perform the actual compilations. If required
  libraries are missing, ``CMake`` will report a corresponding error
  message. In this case, the missing libraries or compilers need to be
  installed before the configuration can be completed.
| To run the ``CMake`` configuration, first create the ``build`` folder
  inside the ``FastChem`` source code folder and switch to the folder:

.. code:: bash

   mkdir build
   cd build

For a basic installation, within the folder run ``CMake``:

.. code:: bash

   cmake ..

If the Python interface should be installed as well, run

.. code:: bash

   cmake -DUSE_PYTHON=ON ..

``CMake`` will also try to locate an ``OpenMP`` library to allow
``FastChem`` to be run in parallel. If it cannot detect the library,
only the single-core version of ``FastChem`` will be compiled. If
``FastChem`` is to be run on MacOS, using ``OpenMP`` might be difficult
since Apple’s ``Clang`` compiler does not directly support ``OpenMP``,
even if the corresponding library has been installed. It might be
possible, though, to install an alternative compiler, for example
``g++``, that supports the use of ``OpenMP``.

After ``CMake`` successfully configured the compilation files,
``FastChem`` can be compiled by running:

.. code:: bash

   make

Upon successful compilation, the executable ``fastchem`` should be
present in the main ``FastChem`` folder. If the optional Python
interface is used, ``FastChem`` will be automatically compiled twice
because the Python version requires different compiler options.

Notes on MacOS
^^^^^^^^^^^^^^

``FastChem`` can be compiled and run on MacOS, but requires some
libraries and apps that are not installed by default. This especially
includes ``CMake``. In order to compile ``FastChem`` on MacOS, the the
prerequisites listed above need to be installed. This can be easily
achieved by, for example, using ``brew``.

In a standard installation of MacOS, no compiler is available. The Apple
version of the Clang compiler can be installed through Xcode and the
command line tools by running

.. code:: bash

   xcode-select --install

in the terminal.

Alternatives (e.g. ``g++``) to the default Clang shipped with MacOS can
also be installed via ``brew``. However, ``CMake`` is not always able to
detect these compilers and will still use Clang. This also applies to
the optional ``OpenMP`` library that allows ``FastChem`` to be run in
parallel. The Clang compiler does not directly support the library, even
if it has been installed via ``brew``.

If the Python interface of ``FastChem`` is used, a corresponding Python
3 installation is also required. By default, MacOS ships only with an
outdated Python 2 version that cannot be used for ``FastChem``. A more
up-to-date version can also be installed by, for example, ``brew``.
However, one has to make sure that the ``python3`` executable and things
like ``pip3`` (to install other required Python modules) actually link
to that version. An alternative way to install and manage different
versions of Python without interference from MacOS’ internal Python
version is ``pyenv``, which can be found under
https://github.com/pyenv/pyenv.

Notes on Windows
^^^^^^^^^^^^^^^^

While in theory ``FastChem`` could be run on Windows if meeting all the
prerequisites, we have never tested the compilation and execution of
``FastChem`` on such a system. In principle, this should be possible
under a virtual Linux environment, such as ``cygwin``, or with the
Windows Subsystem for Linux (WSL) shipped with the newer versions of
Windows 10. However, due to the lack of a Windows system, we are unable
to test this and, therefore, officially at least we cannot support
``FastChem`` running on Windows.

.. _sec:install_python:

Installation of pyFastChem with Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When setting up ``pyFastChem`` with ``PyPI``, it is installed via
``pip``:

.. code:: bash

   pip install pyfastchem

Depending on the Python installation, ``pip`` might need to be replaced
by ``pip3`` in case ``pip`` is linked to Python 2.x.

| This command will download and compile the ``pyFastChem`` package and
  resolve potential dependencies. It is important to note, though, that
  one still has to download the chemistry input data and other Python
  scripts from the ``FastChem`` repository in order to use the package
  properly.
| As an alternative, ``pyFastChem`` can also be directly installed from
  source via the ``setup.py`` script located in the root directory of
  ``FastChem``. The setup is started by

.. code:: bash

   python setup.py install

| assuming that ``python`` points to your Python 3.x executable.
  Otherwise, replace ``python`` with ``python3``. As discussed above,
  using ``setup.py`` will only create the Python module of ``FastChem``,
  not the stand-alone executable.
| In both cases, the ``pyFastChem`` module itself will be installed in
  your local Python package library path and, thus, be available
  throughout your system like any other normal Python package. The
  module’s location and additional module information can be obtained
  via

.. code:: bash

   pip show pyfastchem

The setup script will also try to detect the presence of compiler
support for ``OpenMP`` to run ``FastChem`` calculations in parallel.
This is currently likely to fail in case of MacOS since Apple’s Clang
compiler officially does not support this library. We might adapt the
``setup.py`` script in the future to allow for alternative compilers
under MacOS. 
