# python_coil_calculator
The python code is almost a complete one to one conversion of the code from [Coil64](https://coil32.net/).

Errors that exist in Coil64 will therefore also be present in the code. Furthermore, there are no security mechanisms that check input. To be sure, the required coil calculation in Coil64 should also be tested and compared with the result of the Python code.

USE OF THE CODE AT YOUR OWN RISK!

The Python code is used for easy implementation of automated calculations, which were not possible because of missing API to Coil64.

The Code based on: https://github.com/radioacoustick/Coil64/commit/6532d3ee6c6e5264e533abf201ee716ee4b7c1fd

## usage of the code (100% Python implementation)
Examples for using the code can be found in [example_calc.py](lib_python/example_calc.py) and [example_calc_L.py](lib_python/example_calc_L.py). 

## usage of c-code lib
python setup.py build
pip install .