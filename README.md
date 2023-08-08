# pythonCoil64
The python code is almost a complete one to one conversion of the code from [Coil64](https://coil32.net/).

Errors that exist in Coil64 will therefore also be present in the code. Furthermore, there are no security mechanisms that check input. To be sure, the required coil calculation in Coil64 should also be tested and compared with the result of the Python code.

USE OF THE CODE AT YOUR OWN RISK!

The Python code is used for easy implementation of automated calculations, which were not possible because of missing API to Coil64.

The Code based on: https://github.com/radioacoustick/Coil64/commit/df49e50f94745d08d751946eeb7d2bc01072d62c

## usage of c++-code lib
````
python setup.py build
pip install .
````

For windows users it may be necessary to install a c++ compiler. The library already exists in precompiled form, but in case of errors the library must be completely recompiled.

## Example with images
Tested examples ([example.ipynb](lib_python/example.ipynb)) are available for all functions contained in the library. 

## usage of the code (100% Python implementation)
Examples for using the code can be found in [example_calc.py](lib_python/example_calc.py) and [example_calc_L.py](lib_python/example_calc_L.py). Unfortunately, there are significantly more errors than in the pure C++ version. Therefore please use especially carefully.