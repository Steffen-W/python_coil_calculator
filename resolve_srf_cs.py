# /* resolve_srf_cs.cpp - source text to Coil64 - Radio frequency inductor and choke calculator
# Copyright (C) 2019 Kustarev V.

# This program is free software you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses
# */


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /// SELF-CAPACITANCE & SELF-RESONANCE OF THE ONE-LAYER COIL
# /// http://www.g3ynh.info/zdocs/magnetics/appendix/self_res/self-res.pdf
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

import math
from bessel import *

M_PI = math.pi

e0 = 8.854187818
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def R0tan(x: float) -> float:
    # Calculates R0 Tan(psi) for the sheet helix. D W Knight. Version 1.00, 2016-04-04
    I0 = bessi0(x)
    K0 = bessk0(x)
    I1 = bessi1(x)
    K1 = bessk1(x)
    Result = 59.9584916 * math.sqrt(2 * I0 * K0 * I1 * K1)
    return Result

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def OLLENDF(x: float) -> float:
    # Calculates Ollendorff's function. D W Knight. Version 1.00, 2016-03-16
    I0 = bessi0(x)
    K0 = bessk0(x)
    I1 = bessi1(x)
    K1 = bessk1(x)
    Result = math.sqrt(I0 * K0 / (I1 * K1))
    return Result

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def W82W(x: float) -> float:
    # calculates Nagaoka's coeff. using Wheeler's 1982 eqn (7) as modified by Bob Weaver88.
    # Max error is +/- 21ppM. x = Diam/length. D W Knight, v1.00, June 2012.

    if (x == 0):
        Result = 1
    else:
        zk = 2 / (M_PI * x)
        K0 = 1 / (math.log(8 / M_PI) - 0.5)
        k2 = 24 / (3 * M_PI * M_PI - 16)
        w = -0.47 / pow((0.755 + x), 1.44)
        p = K0 + 3.437 / x + k2 / (x * x) + w
        Result = zk * (math.log(1 + 1 / zk) + 1 / p)

    return Result

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def VFnom(lod: float, ei: float, ex: float) -> float:
    # Calculates nominal helical velocity factor for a free coil at its first SRF.
    # D W Knight, v1.00, 2016-04-14
    # Calls functions W82W(), Ollendf(), R0tan()
    diff = 1
    y = 1

    x = M_PI / (2 * lod)
    x = x + 0.117 * math.exp(-22 * math.exp(-0.75 * lod))
    kL = W82W(1 / lod)
    erad = (ex / 2) * (1 + kL + (ei / ex) * (1 - kL))
    vf0 = 1 / math.sqrt(erad)
    vfh = vf0 * OLLENDF(x)
    Cff = erad * 2 * e0 / (1 + math.log(1 + lod))
    arg = 1.015 + lod * lod
    Caf = e0 * (ei + ex) / math.log(arg + math.sqrt(arg * arg - 1))
    Rotn = vf0 * R0tan(x)
    b = 1000000 * lod / (299.792458 * M_PI * (Cff + Caf) * Rotn)
    a = M_PI / (2 * vfh)
    # solve for VFnom
    z = vfh
    n = 0
    while ((math.fabs(diff) > 1E-9) or (n < 255)):
        y = (1 / a) * math.atan(b / z)
        diff = z - y
        der = -(b / a) / (z * z + b * b)
        deltaz = diff / (der - 1)
        z = z + deltaz
        n += 1

    return y

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def CTDW(ff: float, ei: float, ex: float) -> float:
    # Calculates solenoid time delay capacitance, CT/D [pF/m]. ff is solenoid length / Diameter
    # Quick version using W82W for Nagaoka's coeff. D W Knight. v1.00, 2016-03-16
    kL = W82W(1 / ff)
    kct = 1 / kL - 1
    return 11.27350207 * ex * ff * (1 + kct * (1 + ei / ex) / 2)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def CIAE(ff: float, ei: float, ex: float) -> float:
    # Calculates induced axial E-field component of solenoid self-C, CT/D [pF/m].
    # D W Knight. v1.00, 2016-03-16 . ff is solenoid length / Diameter
    return 17.70837564 * (ei + ex)/math.log(1 + M_PI * M_PI * ff)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /// PUBLIC FUNCTIONS
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def findSRF(lk: float, Dk: float, lw: float) -> float:
    # lk - winding length
    # Dk - winding diameter
    # lw - length of wire
    Vhx = VFnom(lk / Dk, 1, 1)
    Result = Vhx * 299.792458 / (2 * lw)
    return Result

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def find_Cs(p: float, Dk: float, lk: float) -> float:
    # The self-resonance and self-capacitance of solenoid coils D W Knight Version 2 1.00, 4 th May 2016
    # Calculate solenoid self-capacitance, Cl-TDE (9.14) p68
    p = p / 1000
    Dk = Dk / 1000
    lk = lk / 1000
    sinpsi = p / (M_PI * Dk)
    cospsi = math.sqrt(1 - (sinpsi * sinpsi))
    result = CTDW(lk/Dk, 1, 1)/pow(cospsi, 2) + CIAE(lk/Dk, 1, 1)
    return result * Dk
