import matplotlib.pyplot as plt
from inspect import getmembers, isfunction
import numpy as np

import pythonCoil64

pythonCoil64.calc_Onelayer_p()

# xx = np.arange(100)*0.01
# out = []

# for x in xx:
#     temp = pythonCoil64.calc_L_PCB_coil_Square(f=1, t=x)

#     out.append(temp["Q"])

# plt.plot(xx, out)
# plt.show()
