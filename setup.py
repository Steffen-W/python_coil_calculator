from setuptools import setup, Extension
import os

# python setup.py build # to build
# python setup.py bdist # create such a built distribution
# pip install . # to install lib

lib_c_dir = os.path.join(os.getcwd(), "lib_c")

# compile C files:
# g++ -c *.cpp
# ar rvs libcoil64.a *.o

onlyCfiles = ["bessel.cpp", "resolve_q.cpp",
              "resolves.cpp", "resolve_srf_cs.cpp", "coil64_lib.cpp"]

onlyCfiles = [os.path.join(lib_c_dir, f) for f in onlyCfiles]
# wrapper = os.path.join(lib_c_dir, "coil64_lib.cpp")

module = Extension('coil64_lib', sources=onlyCfiles)

setup(
    name='pythonCoil64',
    version='1.0.0',
    author='Steffen-W',
    url='https://github.com/Steffen-W/python_coil_calculator',
    description='Coil64 python API',
    ext_modules=[module],
    packages=['pythonCoil64'],
)
