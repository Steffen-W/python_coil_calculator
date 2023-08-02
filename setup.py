from setuptools import setup, Extension
import os

lib_c_dir = os.path.join(os.getcwd(), "lib_c")
onlyCfiles = [os.path.join(lib_c_dir, f)
              for f in os.listdir(lib_c_dir) if f.endswith(".cpp")]

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
