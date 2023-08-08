from setuptools import setup, Extension
import os
import platform

# to install lib:
# python setup.py build
# pip install .

lib_c_dir = os.path.join(os.getcwd(), "lib_c")
lib_bin_dir = os.path.join(os.getcwd(), "bin")

# coil64_lib.c currently is not working
files = [os.path.join(lib_c_dir, "coil64_lib.cpp")]
if (platform.system() == 'Linux') and (os.path.exists(os.path.join(lib_bin_dir, 'Linux', "libcppcoil64.a"))):       # if static library exist
    module = Extension('coil64_lib', sources=files,
                       library_dirs=[os.path.join('bin', 'Linux')], libraries=[':libcppcoil64.a'])
elif (platform.system() == 'Windows') and (os.path.exists(os.path.join(lib_bin_dir, 'Windows', "libcppcoil64.lib"))):  # if static library exist
    module = Extension('coil64_lib', sources=files,
                       library_dirs=[os.path.join('bin', 'Windows')], libraries=[':libcppcoil64.lib'])
else:  # no precompiled c++ lib exist
    print("no precompiled c++ lib exist!")
    files = ["bessel.cpp", "resolve_q.cpp",
             "resolves.cpp", "resolve_srf_cs.cpp", "coil64_lib.cpp"]
    files = [os.path.join(lib_c_dir, f) for f in files]
    module = Extension('coil64_lib', sources=files)


setup(
    name='pythonCoil64',
    version='1.0.0',
    author='Steffen-W',
    url='https://github.com/Steffen-W/python_coil_calculator',
    description='Coil64 python API',
    ext_modules=[module],
    packages=['pythonCoil64'],
)
