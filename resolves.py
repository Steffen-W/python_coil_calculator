# /* resolves.cpp - source text to Coil64 - Radio frequency inductor and choke calculator
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

import math
from enum import Enum
import numpy as np


class _Elliptic:
    def __init__(self):
        self.Fk = 0.0
        self.Ek = 0.0


class _MagCoreConst:
    def __init__(self):
        self.C1 = 0.0
        self.C2 = 0.0


class _CoilResult:
    def __init__(self):
        self.N = 0.0
        self.sec = 0.0
        self.thd = 0.0
        self.fourth = 0.0
        self.five = 0.0
        self.six = 0
        self.seven = 0.0

    def __str__(self):
        return "N {}, sec {}, thd {}, fourth {}, five {}, six {}, seven {}".format(self.N, self.sec, self.thd, self.fourth, self.five, self.six, self.seven)


Material = {"Al": 0,
            "Cu": 1,
            "Ag": 2,
            "Ti": 3}

M_SQRT2 = math.sqrt(2)
M_PI = math.pi
mu0 = 4e-7 * M_PI

Rho = 0
Chi = 1
Alpha = 2
Dencity = 3
mtrl = [[2.824e-8, 2.21e-5, 0.0039, 2.69808], [1.7241e-8, - 9.56e-6, 0.00393, 8.96],
        [1.59e-8, - 2.63e-5, 0.0038, 10.5], [1.15e-7, 2.4e-6, 0.0042, 7.29]]


def odCalc(id: float) -> float:
    # Calculating the outer diameter (od) of the wire with insulation from the internal diameter (id) without insulation
    M = 0.96344
    b = -0.19861
    od1 = math.exp(M * math.log(id) + b)
    od2 = 1.09 * id
    od = max(od1, od2)
    return od


def find_Archimedean_spiral_length(n: int, a: float) -> float:
    # function to calculate the Archimedean spiral length
    phi = 2 * n * M_PI
    l = (a / (4 * M_PI)) * (phi * math.sqrt(1 + phi * phi) +
                            math.log(phi + math.sqrt(1 + phi * phi)))
    return l


def find_actual_spiral_length(N: int, Din: float, k: float) -> float:
    # Find the spiral length
    ni = math.ceil(Din / (2 * k))
    Lin = find_Archimedean_spiral_length(ni, k)
    n = ni + N
    Lt = find_Archimedean_spiral_length(n, k)
    return (Lt - Lin)


def rosaKm(n: float) -> float:
    # Rosa's round wire mutual inductance correction
    n2 = n * n
    n3 = n2 * n
    n5 = n3 * n2
    n7 = n5 * n2
    n9 = n7 * n2
    return (math.log(2 * M_PI) - 1.5 - math.log(n) / (6 * n) - 0.33084236 / n - 1 / (120 * n3) + 1 / (504 * n5) - 0.0011923 / n7 + 0.0005068 / n9)


def rosaKs(x: float) -> float:
    # Rosa's round wire self inductance correction
    return (1.25 - math.log(2 * x))


def EF(c: float):
    # the complete elliptic integrals of the first and second kind
    a = 1
    b = math.sqrt(1 - math.pow(c, 2))
    E = 1 - math.pow(c, 2) / 2
    i = 1
    while (math.fabs(a - b) > 1e-15):
        a1 = (a + b) / 2
        b1 = math.sqrt(a * b)
        E = E - i * math.pow((a - b) / 2, 2)
        i = 2 * i
        a = a1
        b = b1

    Fk = M_PI / (2 * a)
    Ek = E * Fk
    result = _Elliptic()
    result.Ek = Ek
    result.Fk = Fk
    return result


def Mut(r1: float, r2: float, x: float, g: float) -> float:
    # Mutual inductance of two coaxial circular filaments
    # r1,r2 - radii of the two circular filaments
    # x - distance between the centres of the circular filaments
    # g - Geometric Mean Distance
    l = math.sqrt(math.pow(r2 - r1, 2) + math.pow(x, 2))
    c = 2 * math.sqrt(r1 * r2) / \
        math.sqrt(math.pow(r1 + r2, 2) + math.pow(l - g, 2))
    Ec = EF(c)
    result = -0.004 * M_PI * \
        math.sqrt(r1 * r2) * ((c - 2 / c) * Ec.Fk + (2 / c) * Ec.Ek)
    return result


def SelfInductanceStraightWire(l: float, dw: float) -> float:
    r = 0.5 * dw * math.exp(0.25)
    result = 0.002 * (l * math.log((l + math.sqrt(l * l + r * r)
                                    ) / r) - math.sqrt(l * l + r * r) + l / 4 + r)
    return result


def MutInductanceStraightWire(L1: float, l2: float, D: float) -> float:
    tmpLength = 0
    L1 = L1 / 2
    l2 = l2 / 2
    if (L1 < l2):
        tmpLength = L1
        L1 = l2
        l2 = tmpLength

    result = 0.002 * (2 * L1 * math.log((L1 + l2 + math.sqrt((L1 + l2) * (L1 + l2) + D * D)) / D) + (L1 - l2) *
                      math.log((L1 + l2 + math.sqrt((L1 + l2) * (L1 + l2) + D * D)) / (L1 - l2 + math.sqrt((L1 - l2) * (L1 - l2) + D * D))) +
                      math.sqrt((L1 - l2) * (L1 - l2) + D * D) - math.sqrt((L1 + l2) * (L1 + l2) + D * D))
    return result


def Ingrnd(phi: float, kphitheta: float, sinpsi: float, cos2psi: float, rr: float, y: float) -> float:
    # by Robert Weaver from http://electronbunker.ca/eb/CalcMethods2d.html
    # Integrand function called by HeliCoilS()
    result = (1 + cos2psi * (math.cos(kphitheta) - 1)) / math.sqrt(2 * rr * (1 - math.cos(kphitheta)) + (sinpsi * phi - y) *
                                                                   (sinpsi * phi - y))
    return result


def HeliCoilS(Lw: float, psi: float, r: float, dw: float, w: float, t: float, isRoundWire: bool, accuracy: int) -> float:
    # by Robert Weaver from http://electronbunker.ca/eb/CalcMethods2d.html (Version 1.0, 2011-03-25)
    # edited by Valery Kustarev 2018-12-16

    # Uses helical filament mutual inductance formula
    # evaluated using Simpson's rule, and conductor gmd
    # Lw = total length of wire
    # psi = pitch angle of winding
    # r = radius of winding
    # dw = wire diameter
    # MaxErr = max allowable error (set to 10000)
    #
    # px = p/(2pi) = lCOIL/(2piN)
    # psi = Arctan(px/(piD)) = Arctan(p/(2pi^2*D))
    # lW = piND/cos psi
    # where:
    # D is coil diameter = 2r
    # p is centre to centre turns spacing
    # lCOIL is length of coil.
    # If Lw>2*pi*r, check that pitch angle >= psi-c (close wound pitch)

    err = math.ceil(accuracy / 2)
    MaxErr = math.pow(10, -err)
    Integral = 0
    #  grandtotal = 0
    if (Lw > 2 * M_PI * r):
        if (isRoundWire):
            sinpsic = dw / (2 * M_PI * r)
        else:
            sinpsic = w / (2 * M_PI * r)

        psic = math.atan(sinpsic / math.sqrt(1 - sinpsic * sinpsic))
        if (psi < psic):
            #     pitch angle is too small,#      so set value of function to an illegal value and exit
            return -1

    # gmd of solid round conductor. Other values may be substituted
    # for different conductor geometries
    if (isRoundWire):
        g = math.exp(-0.25) * dw / 2
    else:
        g = 0.2235*(w + t)

    rr = r * r
    psio = 0.5 * M_PI - psi
    # Calculate Filament 2 offset angle
    # Trap for psi=0 condition in which case ThetaO=0 and Y0=g
    # Trap for psio=0 condition in which case use simplified formula for ThetaO and Y0=0
    # which happens with circular (non-helical) filament
    if (psi == 0):
        ThetaO = 0
        Y0 = g
    elif (psio == 0):
        cosThetaO = 1 - (g * g / (2 * rr))
        ThetaO = - \
            math.fabs(math.atan(math.sqrt(1 - cosThetaO * cosThetaO) / cosThetaO))
        Y0 = 0
    else:
        #  Use Newton-Raphson method
        k1 = (g * g) / (2 * r * r) - 1
        k2 = math.tan(psio)
        k2 = 0.5 * k2 * k2
        t1 = g / r * math.sin(psi)

        t0 = t1
        t1 = t0 - (k1 + math.cos(t0) - k2 * t0 * t0) / \
            (-math.sin(t0) - 2 * k2 * t0)
        while (math.fabs(t1 - t0) > 1e-12):
            t0 = t1
            t1 = t0 - (k1 + math.cos(t0) - k2 * t0 * t0) / \
                (-math.sin(t0) - 2 * k2 * t0)

        ThetaO = -math.fabs(t1)
        #   Calculate Filament 2 Y-offset, using formula (29)
        Y0 = math.sqrt(g * g - 2 * rr * (1 - math.cos(ThetaO)))

    # Psi constants
    c2s = math.cos(psi) * math.cos(psi)
    ss = math.sin(psi)
    k = math.cos(psi) / r
    # Start of Simpson's rule code
    a = 0
    b = Lw / 32768
    if (b > Lw):
        b = Lw

    grandtotal = 0
    while (a < Lw):
        dx = b - a
        m = 1
        CurrentErr = 2 * MaxErr
        kat = k * a
        kbt = k * b
        aaa = (Lw - a) * (Ingrnd(-a, -kat - ThetaO, ss, c2s, rr, Y0))
        bbb = Ingrnd(a, kat - ThetaO, ss, c2s, rr, Y0)
        ccc = (Lw - b) * (Ingrnd(-b, -kbt - ThetaO, ss, c2s, rr, Y0))
        ddd = Ingrnd(b, kbt - ThetaO, ss, c2s, rr, Y0)
        Sum2 = aaa + bbb + ccc + ddd

        #   Initialize LastResult to trapezoidal area for test purposes

        LastIntg = Sum2 / 2 * dx
        while ((CurrentErr > MaxErr) or (m < 1024)):
            m = 2 * m
            dx = dx / 2
            Sum = 0
            max = round(m/2)
            for i in range(1, max + 1):
                phi = 2 * i * dx + a
                kpt = k * phi
                Sum = Sum + (Lw - phi) * (Ingrnd(-phi, -kpt - ThetaO, ss,
                                                 c2s, rr, Y0) + Ingrnd(phi, kpt - ThetaO, ss, c2s, rr, Y0))

            Integral = (4 * (Sum) + Sum2) * dx / 3
            CurrentErr = math.fabs((Integral) / (LastIntg) - 1)
            LastIntg = Integral
            Sum2 = Sum2 + Sum * 2

        grandtotal += Integral
        a = b
        b = b * 2
        if (b > Lw):
            b = Lw

    result = 1e-1 * grandtotal
    return result


def solveHelicalInductance(N: float, _p: float, _Dk: float, _dw: float, _w: float, _t: float, isRoundWire: bool, accuracy: int) -> float:
    dw = 0
    w = 0
    t = 0

    p = _p / 1000
    Dk = _Dk / 1000

    if (isRoundWire):
        dw = _dw / 1000
    else:
        w = _w / 1000
        t = _t / 1000

    # psi= math.atan(p/(2*pi*pi*Dk))
    sinpsi = p / (M_PI * Dk)
    psi = math.atan(sinpsi / math.sqrt(1 - sinpsi * sinpsi))
    lW = M_PI * N * Dk / math.cos(psi)
    if (isRoundWire):
        Result = HeliCoilS(lW, psi, Dk / 2, dw, 0, 0, isRoundWire, accuracy)
    else:
        Result = HeliCoilS(lW, psi, Dk / 2, 0, w, t, isRoundWire, accuracy)

    return Result, lW


def deriveOneLayerPoligonalN(Dk: float, dw: float, p: float, n: float, I: float, accuracy: int) -> float:

    k = 2
    N_min = 0
    rA = math.sqrt((1 / M_PI) * (0.5 * n * math.pow(0.5 * Dk, 2)
                   * math.sin(2 * M_PI / n)))
    rP = (0.5 / M_PI) * (Dk * n * math.sin(M_PI / n))
    Kw = math.sqrt(1 / 369.0)
    iDk = 2 * (((Kw * math.pow(rP, 2)) + ((2 - Kw) * math.pow(rA, 2))) / (2 * rA))
    N = math.sqrt(I / (0.0002 * M_PI * iDk * (math.log(1 + M_PI / (2 * k)) + 1 /
                  (2.3004 + 3.437 * k + 1.763 * k * k - 0.47 / math.pow((0.755 + 1 / k), 1.44)))))
    res = _CoilResult()
    getOneLayerI_Poligonal(Dk, dw, p, N, n, res, accuracy)
    Ind = res.sec
    while (Ind < I):
        N_min = N
        N_max = 2 * N
        N = (N_max + N_min) / 2
        getOneLayerI_Poligonal(Dk, dw, p, N, n, res, accuracy)
        Ind = res.sec

    N_max = N
    while (math.fabs(1 - (Ind / I)) > 0.001):
        N = (N_min + N_max) / 2
        getOneLayerI_Poligonal(Dk, dw, p, N, n, res, accuracy)
        Ind = res.sec
        if (Ind > I):
            N_max = N
        else:
            N_min = N

    lw = res.thd
    return N, lw, iDk


def getOneLayerN_Poligonal(I: float, Dk: float, dw: float,  p: float, n: float, result: _CoilResult, accuracy: int) -> float:

    N, lw, iDk = deriveOneLayerPoligonalN(Dk, dw, p, n, I, accuracy)
    result.sec = p * N
    result.thd = lw
    result.seven = iDk
    return N


def getOneLayerI_Poligonal(Dk: float, dw: float, p: float, N: float, n: float, result: _CoilResult, accuracy: int):

    lk = N * p
    rA = math.sqrt((1 / M_PI) * (0.5 * n * math.pow(0.5 * Dk, 2)
                   * math.sin(2.0 * M_PI / n)))
    rP = (0.5 / M_PI) * (Dk * n * math.sin(M_PI / n))
    Kw = math.sqrt(1 / (1 + 368.0 * (lk / Dk)))
    iDk = 2.0 * (((Kw * math.pow(rP, 2)) +
                 ((2.0 - Kw) * math.pow(rA, 2))) / (2.0 * rA))
    result.sec, lw = solveHelicalInductance(
        N, p, iDk, dw, 0, 0, True, accuracy)
    result.thd = lw
    result.seven = iDk


def getFerriteCoreMagConst(l1: float, l2: float, l3: float, l4: float, l5: float,
                           A1: float, A2: float, A3: float, A4: float, A5: float):
    # auxiliary function to get the constants C1 & C2 of a ferrite core with the close magnetic circuit
    sum11 = l1 / A1
    sum21 = l1 / (A1 * A1)
    sum12 = l2 / A2
    sum22 = l2 / (A2 * A2)
    sum13 = l3 / A3
    sum23 = l3 / (A3 * A3)
    sum14 = l4 / A4
    sum24 = l4 / (A4 * A4)
    sum15 = l5 / A5
    sum25 = l5 / (A5 * A5)
    result = _MagCoreConst()
    result.C1 = sum11 + sum12 + sum13 + sum14 + sum15
    result.C2 = sum21 + sum22 + sum23 + sum24 + sum25
    return result

# PUBLIC FUNCTIONS REALIZATION


def getOneLayerN_withRoundWire(Dk: float, dw: float, p: float, I: float, accuracy: int) -> float:

    k = 2
    N_min = 0

    N = math.sqrt(I / (0.0002 * M_PI * Dk * (math.log(1 + M_PI / (2 * k)) + 1 /
                  (2.3004 + 3.437 * k + 1.763 * k * k - 0.47 / math.pow((0.755 + 1 / k), 1.44)))))
    ind, lw = solveHelicalInductance(N, p, Dk, dw, 0, 0, True, accuracy)
    while (ind < I):
        N_min = N
        N_max = 2 * N
        N = (N_max + N_min) / 2
        ind, lw = solveHelicalInductance(N, p, Dk, dw, 0, 0, True, accuracy)

    N_max = N
    while (math.fabs(1 - (ind / I)) > 0.001):
        N = (N_min + N_max) / 2
        ind, lw = solveHelicalInductance(N, p, Dk, dw, 0, 0, True, accuracy)
        if (ind > I):
            N_max = N
        else:
            N_min = N

    return N, lw


def getOneLayerN_byWindingLength(D, L: float, I: float, result: _CoilResult, accuracy: int) -> float:
    dw = 0
    lTmp = 0
    N = 0
    dw_max = 0.25 * D
    dw_min = 0
    i = 0
    while (abs(1 - lTmp/L) > 0.05):
        dw = (dw_min + dw_max) / 2
        k = odCalc(dw)
        Dk = D + k
        N, lw = getOneLayerN_withRoundWire(Dk, dw, k, I, accuracy)
        lTmp = N * k + k
        if (lTmp > L):
            dw_max = dw
        else:
            dw_min = dw

        i += 1
        if (i > 500):
            return 0

    result.sec = lw
    result.five = dw
    return N


def getOneLayerI_withRoundWire(Dk: float, dw: float, p: float, N: float, accuracy: int) -> float:
    ret, lw = solveHelicalInductance(N, p, Dk, dw, 0, 0, True, accuracy)
    return ret, lw


def getOneLayerN_withRectWire(Dk: float, w: float, t: float, p: float, I: float, accuracy: int) -> float:

    k = 2
    N_min = 0

    N = math.sqrt(I / (0.0002 * M_PI * Dk * (math.log(1 + M_PI / (2 * k)) + 1 /
                  (2.3004 + 3.437 * k + 1.763 * k * k - 0.47 / math.pow((0.755 + 1 / k), 1.44)))))
    ind, lw = solveHelicalInductance(N, p, Dk, 0, w, t, False, accuracy)
    while (ind < I):
        N_min = N
        N_max = 2 * N
        N = (N_max + N_min) / 2
        ind, lw = solveHelicalInductance(N, p, Dk, 0, w, t, False, accuracy)

    N_max = N
    while (math.fabs(1 - (ind / I)) > 0.001):
        N = (N_min + N_max) / 2
        ind, lw = solveHelicalInductance(N, p, Dk, 0, w, t, False, accuracy)
        if (ind > I):
            N_max = N
        else:
            N_min = N

    return N, lw


def getOneLayerI_withRectWire(Dk: float, w: float, t: float, p: float, N: float, accuracy: int) -> float:
    res, lw = solveHelicalInductance(N, p, Dk, 0, w, t, False, accuracy)
    return res, lw


def getMultiLayerN(I: float, D: float, dw: float, k: float, lk: float, gap: float, Ng: int, result: _CoilResult):
    n_g = 0
    jg = 0
    D = D / 10
    lk = lk / 10
    dw = dw / 10
    k = k / 10
    gap = gap / 10
    if (Ng == -1):
        gap = 0

    Ltotal = 0  # initialize variable of total self-inductance
    nLayer = 1
    lw = 0
    r0 = (D + k) / 2
    N = 0
    Nl = math.floor(lk / k)  # number of turns in layer
    g = math.exp(-0.25) * dw / 2
    while (Ltotal < I):  # start calculation loop increasing N-turns to reach requiring inductance (I)

        N += 1
        Nc = (N - 1) % Nl  # position of N-turn in layer
        nLayer = math.floor((N - 1) / Nl)  # current layer for N-turn
        if (((nLayer % Ng) == 0) and (nLayer > 0)):
            n_g = gap
        else:
            n_g = 0

        nx = Nc * k  # x-offset of turn
        ny = r0 + k * nLayer + n_g  # y-offset of turn
        Lns = Mut(ny, ny, g, 0)  # self inductance of current turn
        lw = lw + 2 * M_PI * ny  # length of wire with the current turn
        # start calculation loop of the mutual inductance - current turn (N) + all another turns (j)
        M = 0
        if (N > 1):
            for j in range(N, 1, -1):
                Jc = (j - 2) % Nl
                jx = Jc * k
                jLayer = math.floor((j - 2) / Nl)
                if (((jLayer % Ng) == 0) and (jLayer > 0)):
                    jg = gap
                else:
                    jg = 0

                jy = r0 + k * jLayer + jg
                M = M + 2 * Mut(ny, jy, nx - jx, g)  # mutual inductance
                # between current
                # N-turn and j-turn

        # total summary inductance (adding self-inductance and mutual inductance of current N-turn)
        Ltotal += Lns + M

    Resistivity = mtrl[Material["Cu"]][Rho]*1e2
    R = (Resistivity * lw * 4) / (M_PI * dw * dw)  # resistance of the wire
    lw0 = lw / 100
    NLayer = nLayer + 1
    NumberInterLayer = math.floor(nLayer / Ng)
    c = NLayer * k * 10 + NumberInterLayer * gap * 10
    result.N = R
    result.sec = lw0
    result.thd = NLayer
    result.fourth = c
    result.five = NumberInterLayer
    result.six = N
    # return: N(turns) R(resistance) Ohm lw0(length of wire) m
    # NLayer (Number of layer) c (Winding thickness) mm
    # NumberInterLayer (Number of inter-layers)


def getMultiLayerN_rectFormer(Ind: float, a: float, b: float, l: float, dw: float, k: float, result: _CoilResult):
    # Calculation formulas of multilayer inductor with rectangular former https://coil32.net/multilayer-rectangular.html

    a = a / 10
    b = b / 10
    l = l / 10
    dw = dw / 10
    k = k / 10
    n = 0
    Ltotal = 0
    a0 = a + k
    b0 = b + k
    lw = 0
    nLayer = 0
    Nl = math.floor(l / k)  # Number of turns in layer
    while (Ltotal < Ind):
        n += 1
        Nc = (n - 1) % Nl  # Position of the turn on x
        nLayer = math.floor((n - 1) / Nl)  # Position of the turn on y
        nx = Nc * k  # x-offset of current turn
        ny = nLayer * k  # y-offset of current turn
        lengthNa = a0 + 2 * k * (nLayer)
        # lenght of straight conductor of current turn (side a)
        Rdc = (0.0175 * lw * 1E-4 * 4) / (M_PI * dw * dw)
        # lenght of straight conductor of current turn (side b)
        lengthNb = b0 + 2 * k * (nLayer)
        lw += 2 * (a0 + b0 + 2 * k * (nLayer))
        # half of self-inductance of the current turn
        Ladd = SelfInductanceStraightWire(
            lengthNa, dw) + SelfInductanceStraightWire(lengthNb, dw)
        # distance to opposite cunductor of the same turn (side b)
        Db = 2 * ny + a0
        # distance to opposite cunductor of the same turn (side a)
        Da = 2 * ny + b0
        Lsub = 2 * MutInductanceStraightWire(
            lengthNa, lengthNa, Da) + 2 * MutInductanceStraightWire(lengthNb, lengthNb, Db)
        # half mutual inductance with opposite conductor of current turn
        Madd = 0
        Msub = 0
        if (n > 1):
            for j in range(n, 1, -1):
                Jc = (j - 2) % Nl  # position of previous turn on x
                jx = Jc * k  # x-offset of previous turn
                jLayer = math.floor((j - 2) / Nl)  # Position of the turn on y
                jy = k * jLayer  # y-offset of previous turn
                # lenght of straight conductor of previous turn (side a)
                lengthJa = a0 + 2 * k * (nLayer + 1)
                # lenght of straight conductor of previous turn (side b)
                lengthJb = b0 + 2 * k * (nLayer + 1)
                # distance to in-phase straight conductor of previous turn
                D = math.sqrt(math.pow(nx - jx, 2) + math.pow(ny - jy, 2))
                Madd = Madd + 2 * MutInductanceStraightWire(
                    lengthNa, lengthJa, D) + 2 * MutInductanceStraightWire(lengthNb, lengthJb, D)
                # half mutual inductance with in-phase conductor in previous turn
                Db = math.sqrt(math.pow(nx - jx, 2) +
                               math.pow(ny + jy + a0, 2))
                # distance to opposite cunductor between the current turn and previous (side b)
                Da = math.sqrt(math.pow(nx - jx, 2) +
                               math.pow(ny + jy + b0, 2))
                # distance to opposite cunductor between the current turn and previous (side a)
                Msub = Msub + 2 * MutInductanceStraightWire(
                    lengthNa, lengthJa, Da) + 2 * MutInductanceStraightWire(lengthNb, lengthJb, Db)
                # half mutual inductance with opposite conductor in previous turn

        Ltotal += 2 * (Ladd - Lsub + Madd - Msub)
        Ks = rosaKs(k / dw)
        Km = rosaKm(n)
        Lcor = 0.0002 * M_PI * (a + b) * n * (Ks + Km)
        Ltotal -= Lcor

    Resistivity = mtrl[Material["Cu"]][Rho]*1e2
    Rdc = (Resistivity * lw * 4) / (M_PI * dw * dw)
    result.N = n  # number of turns
    result.sec = nLayer + 1  # number of layers
    result.thd = lw * 0.01  # length of wire
    result.fourth = Rdc  # resistance to DC
    result.five = (nLayer + 1) * k * 10  # coil thickness


def getMultiLayerI_rectFormer_byN(N: float, a: float, b: float, l: float, dw: float, k: float, result: _CoilResult):

    a = a / 10
    b = b / 10
    l = l / 10
    dw = dw / 10
    k = k / 10
    Ltotal = 0
    a0 = a + k
    b0 = b + k
    lw = 0
    nLayer = 0
    Nl = math.floor(l / k)  # Number of turns in layer
    for n in range(1, N + 1):
        Nc = (n - 1) % Nl  # Position of the turn on x
        nLayer = math.floor((n - 1) / Nl)  # Position of the turn on y
        nx = Nc * k  # x-offset of current turn
        ny = nLayer * k  # y-offset of current turn
        lengthNa = a0 + 2 * k * (nLayer)
        # lenght of straight conductor of current turn (side b)
        lengthNb = b0 + 2 * k * (nLayer)
        lw += 2 * (a0 + b0 + 2 * k * (nLayer))
        # half of self-inductance of the current turn
        Ladd = SelfInductanceStraightWire(
            lengthNa, dw) + SelfInductanceStraightWire(lengthNb, dw)
        # distance to opposite cunductor of the same turn (side b)
        Db = 2 * ny + a0
        # distance to opposite cunductor of the same turn (side a)
        Da = 2 * ny + b0
        Lsub = 2 * MutInductanceStraightWire(
            lengthNa, lengthNa, Da) + 2 * MutInductanceStraightWire(lengthNb, lengthNb, Db)
        # half mutual inductance with opposite conductor of current turn
        Madd = 0
        Msub = 0
        if (n > 1):
            for j in range(n, 1, -1):
                Jc = (j - 2) % Nl  # position of previous turn on x
                jx = Jc * k  # x-offset of previous turn
                jLayer = math.floor((j - 2) / Nl)  # Position of the turn on y
                jy = k * jLayer  # y-offset of previous turn
                # lenght of straight conductor of previous turn (side a)
                lengthJa = a0 + 2 * k * (nLayer + 1)
                # lenght of straight conductor of previous turn (side b)
                lengthJb = b0 + 2 * k * (nLayer + 1)
                # distance to in-phase straight conductor of previous turn
                D = math.sqrt(math.pow(nx - jx, 2) + math.pow(ny - jy, 2))
                Madd = Madd + 2 * MutInductanceStraightWire(
                    lengthNa, lengthJa, D) + 2 * MutInductanceStraightWire(lengthNb, lengthJb, D)
                # half mutual inductance with in-phase conductor in previous turn
                Db = math.sqrt(math.pow(nx - jx, 2) +
                               math.pow(ny + jy + a0, 2))
                # distance to opposite cunductor between the current turn and previous (side b)
                Da = math.sqrt(math.pow(nx - jx, 2) +
                               math.pow(ny + jy + b0, 2))
                # distance to opposite cunductor between the current turn and previous (side a)
                Msub = Msub + 2 * MutInductanceStraightWire(
                    lengthNa, lengthJa, Da) + 2 * MutInductanceStraightWire(lengthNb, lengthJb, Db)
                # half mutual inductance with opposite conductor in previous turn

        Ltotal += 2 * (Ladd - Lsub + Madd - Msub)
        Ks = rosaKs(k / dw)
        Km = rosaKm(n)
        Lcor = 0.0002 * M_PI * (a + b) * n * (Ks + Km)
        Ltotal -= Lcor

    result.N = Ltotal  # inductance
    result.sec = nLayer + 1  # number of layers
    result.thd = lw * 0.01  # length of wire
    result.five = (nLayer + 1) * k * 10  # coil thickness


def getMultiLayerI_byN(D: float, lk: float, dw: float, k: float, N: float, result: _CoilResult):
    D = D / 10
    lk = lk / 10
    dw = dw / 10
    k = k / 10
    nLayer = 0
    Ltotal = 0  # initialize variable of total self-inductance
    lw = 0
    r0 = (D + k) / 2
    Nl = math.floor(lk / k)
    g = math.exp(-0.25) * dw / 2
    for w in range(1, N + 1):
        Nc = (w - 1) % Nl
        nLayer = math.floor((w - 1) / Nl)
        nx = Nc * k
        ny = r0 + k * nLayer
        Lns = Mut(ny, ny, g, 0)
        # self inductance of current turn
        lw = lw + 2 * M_PI * ny
        M = 0
        if (w > 1):
            for j in range(w, 1, -1):
                Jc = (j - 2) % Nl
                jx = Jc * k
                jLayer = math.floor((j - 2) / Nl)
                jy = r0 + k * jLayer
                M = M + 2 * Mut(ny, jy, nx - jx, g)

        Ltotal += Lns + M

    Resistivity = mtrl[Material["Cu"]][Rho]*1e2
    Rdc = (Resistivity * lw * 4) / (M_PI * dw * dw)
    result.N = Ltotal  # inductance value
    result.sec = nLayer + 1  # number of layers
    result.thd = lw * 0.01  # length of wire
    result.fourth = Rdc  # resistance to DC
    result.five = (nLayer + 1) * k * 10  # coil thickness


def getMultiLayerI(D: float, lk: float, dw: float, k: float, c: float, gap: float,  Ng: int, result: _CoilResult):
    n_g = 0
    jg = 0
    ind1 = 0
    D = D / 10
    lk = lk / 10
    c = c / 10
    nTmp = 0
    bTmp = 0
    dw = dw / 10
    k = k / 10
    gap = gap / 10
    if (Ng == -1):
        gap = 0

    Ltotal = 0  # initialize variable of total self-inductance
    lw = 0
    r0 = (D + k) / 2
    n = 0
    Nl = math.floor(lk / k)
    g = math.exp(-0.25) * dw / 2
    while (bTmp < (c + k)):
        n += 1
        Nc = (n - 1) % Nl
        nLayer = math.floor((n - 1) / Nl)
        if (((nLayer % Ng) == 0) and (nLayer > 0)):
            n_g = gap
        else:
            n_g = 0

        nx = Nc * k
        ny = r0 + k * nLayer + n_g
        Lns = Mut(ny, ny, g, 0)
        # self inductance of current turn
        lw = lw + 2 * M_PI * ny
        M = 0
        if (n > 1):
            for j in range(n, 1, -1):
                Jc = (j - 2) % Nl
                jx = Jc * k
                jLayer = math.floor((j - 2) / Nl)
                if (((jLayer % Ng) == 0) and (jLayer > 0)):
                    jg = gap
                else:
                    jg = 0

                jy = r0 + k * jLayer + jg
                M = M + 2 * Mut(ny, jy, nx - jx, g)

        Ltotal += Lns + M
        if (nTmp < c):
            nTmp = (nLayer + 1) * k
            ind1 = Ltotal

        bTmp = (nLayer + 1) * k

    N1 = n - Nl
    N2 = n + Nl
    ind2 = Ltotal
    result.N = ind1
    result.sec = ind2
    result.thd = N1
    result.fourth = N2


def getMultiLayerI_rectFormer(a: float, b: float, l: float, c: float, dw: float, k: float, result: _CoilResult):

    cTmp = 0
    nTmp = 0
    ind1 = 0

    a = a / 10
    b = b / 10
    l = l / 10
    c = c / 10
    dw = dw / 10
    k = k / 10
    n = 0
    Ltotal = 0
    a0 = a + k
    b0 = b + k
    lw = 0
    nLayer = 0
    Nl = math.floor(l / k)  # Number of turns in layer
    while (cTmp < (c + k)):
        n += 1
        Nc = (n - 1) % Nl  # Position of the turn on x
        nLayer = math.floor((n - 1) / Nl)  # Position of the turn on y
        nx = Nc * k  # x-offset of current turn
        ny = nLayer * k  # y-offset of current turn
        # lenght of straight conductor of current turn (side a)
        lengthNa = a0 + 2 * k * (nLayer)
        # lenght of straight conductor of current turn (side b)
        lengthNb = b0 + 2 * k * (nLayer)
        lw = lw + 2 * (a0 + b0 + 2 * k * (nLayer))
        # half of self-inductance of the current turn
        Ladd = SelfInductanceStraightWire(
            lengthNa, dw) + SelfInductanceStraightWire(lengthNb, dw)
        # distance to opposite cunductor of the same turn (side b)
        Db = 2 * ny + a0
        # distance to opposite cunductor of the same turn (side a)
        Da = 2 * ny + b0
        Lsub = 2 * MutInductanceStraightWire(
            lengthNa, lengthNa, Da) + 2 * MutInductanceStraightWire(lengthNb, lengthNb, Db)
        # half mutual inductance with opposite conductor of current turn
        Madd = 0
        Msub = 0
        if (n > 1):
            for j in range(n, 1, -1):
                Jc = (j - 2) % Nl  # position of previous turn on x
                jx = Jc * k  # x-offset of previous turn
                jLayer = math.floor((j - 2) / Nl)  # Position of the turn on y
                jy = k * jLayer  # y-offset of previous turn
                # lenght of straight conductor of previous turn (side a)
                lengthJa = a0 + 2 * k * (nLayer + 1)
                # lenght of straight conductor of previous turn (side b)
                lengthJb = b0 + 2 * k * (nLayer + 1)
                # distance to in-phase straight conductor of previous turn
                D = math.sqrt(math.pow(nx - jx, 2) + math.pow(ny - jy, 2))
                Madd = Madd + 2 * MutInductanceStraightWire(
                    lengthNa, lengthJa, D) + 2 * MutInductanceStraightWire(lengthNb, lengthJb, D)
                # half mutual inductance with in-phase conductor in previous turn
                Db = math.sqrt(math.pow(nx - jx, 2) +
                               math.pow(ny + jy + a0, 2))
                # distance to opposite cunductor between the current turn and previous (side b)
                Da = math.sqrt(math.pow(nx - jx, 2) +
                               math.pow(ny + jy + b0, 2))
                # distance to opposite cunductor between the current turn and previous (side a)
                Msub = Msub + 2 * MutInductanceStraightWire(
                    lengthNa, lengthJa, Da) + 2 * MutInductanceStraightWire(lengthNb, lengthJb, Db)
                # half mutual inductance with opposite conductor in previous turn

        Ltotal += 2 * (Ladd - Lsub + Madd - Msub)
        Ks = rosaKs(k / dw)
        Km = rosaKm(n)
        Lcor = 0.0002 * M_PI * (a + b) * n * (Ks + Km)
        Ltotal -= Lcor
        if (nTmp < c):
            nTmp = (nLayer + 1) * k
            ind1 = Ltotal

        cTmp = (nLayer + 1) * k

    N1 = n - Nl
    N2 = n + Nl
    ind2 = Ltotal
    result.N = ind1
    result.sec = ind2
    result.thd = N1
    result.fourth = N2


def getMultiLayerI_fromResistance(D: float, lk: float, c: float, k: float, Rm: float, result: _CoilResult):

    nLayer = 0
    aWire = [[0.06, 0.075, 0.09, 0.085, 0.09],
             [0.063, 0.078, 0.09, 0.085, 0.09], [0.07, 0.084, 0.092,
                                                 0.092, 0.1], [0.071, 0.088, 0.095, 0.095, 0.1],
             [0.08, 0.095, 0.105, 0.105, 0.11], [0.09, 0.105, 0.12,
                                                 0.115, 0.12], [0.1, 0.122, 0.13, 0.125, 0.13],
             [0.112, 0.134, 0.14, 0.125, 0.14], [0.12, 0.144, 0.15,
                                                 0.145, 0.15], [0.125, 0.149, 0.155, 0.15, 0.155],
             [0.13, 0.155, 0.16, 0.155, 0.16], [0.14, 0.165, 0.17,
                                                0.165, 0.17], [0.15, 0.176, 0.19, 0.18, 0.19],
             [0.16, 0.187, 0.2, 0.19, 0.2], [0.17, 0.197, 0.21,
                                             0.2, 0.21], [0.18, 0.21, 0.22, 0.21, 0.22],
             [0.19, 0.22, 0.23, 0.22, 0.23], [0.2, 0.23, 0.24,
                                              0.23, 0.24], [0.21, 0.24, 0.25, 0.25, 0.25],
             [0.224, 0.256, 0.27, 0.26, 0.27], [0.236, 0.26, 0.285,
                                                0.27, 0.28], [0.25, 0.284, 0.3, 0.275, 0.3],
             [0.265, 0.305, 0.315, 0.305, 0.31], [0.28, 0.315,
                                                  0.33, 0.315, 0.33], [0.3, 0.34, 0.35, 0.34, 0.34],
             [0.315, 0.35, 0.365, 0.352, 0.36], [0.335, 0.375, 0.385,
                                                 0.375, 0.38], [0.355, 0.395, 0.414, 0.395, 0.41],
             [0.38, 0.42, 0.44, 0.42, 0.44], [0.4, 0.44, 0.46,
                                              0.442, 0.46], [0.425, 0.465, 0.485, 0.47, 0.47],
             [0.45, 0.49, 0.51, 0.495, 0.5], [0.475, 0.525, 0.545,
                                              0.495, 0.53], [0.5, 0.55, 0.57, 0.55, 0.55],
             [0.53, 0.58, 0.6, 0.578, 0.6], [0.56, 0.61, 0.63,
                                             0.61, 0.62], [0.6, 0.65, 0.67, 0.65, 0.66],
             [0.63, 0.68, 0.7, 0.68, 0.69], [0.67, 0.72, 0.75,
                                             0.72, 0.75], [0.71, 0.76, 0.79, 0.77, 0.78],
             [0.75, 0.81, 0.84, 0.81, 0.83], [0.8, 0.86, 0.89,
                                              0.86, 0.89], [0.85, 0.91, 0.94, 0.91, 0.94],
             [0.9, 0.96, 0.99, 0.96, 0.99], [0.93, 0.99, 1.02,
                                             0.99, 1.02], [0.95, 1.01, 1.04, 1.02, 1.04],
             [1.0, 1.07, 1.1, 1.07, 1.11], [1.06, 1.13, 1.16,
                                            1.14, 1.16], [1.08, 1.16, 1.19, 1.16, 1.19],
             [1.12, 1.19, 1.22, 1.2, 1.23], [1.18, 1.26, 1.28,
                                             1.26, 1.26], [1.25, 1.33, 1.35, 1.33, 1.36],
             [1.32, 1.4, 1.42, 1.4, 1.42], [1.4, 1.48, 1.51,
                                            1.48, 1.51], [1.45, 1.53, 1.56, 1.53, 1.56],
             [1.5, 1.58, 1.61, 1.58, 1.61], [1.56, 1.63, 1.67,
                                             1.64, 1.67], [1.6, 1.68, 1.71, 1.68, 1.71],
             [1.7, 1.78, 1.81, 1.78, 1.81], [1.74, 1.82, 1.85,
                                             1.82, 1.85], [1.8, 1.89, 1.92, 1.89, 1.92],
             [1.9, 1.99, 2.02, 1.99, 2.02], [2.0, 2.1, 2.12,
                                             2.1, 2.12], [2.12, 2.21, 2.24, 2.22, 2.24],
             [2.24, 2.34, 2.46, 2.34, 2.46], [2.36, 2.46, 2.48, 2.36, 2.48], [2.5, 2.6, 2.63, 2.6, 2.62]]

    D = D / 10
    lk = lk / 10
    k = k / 10
    c = c / 10
    bTmp = 0
    nTmp = 0
    for z in range(67):
        dw = aWire[z][0] / 10
        Ltotal = 0
        # initialize variable of total self-inductance
        lw = 0
        r0 = (D + k) / 2
        n = 0
        Nl = math.floor(lk / k)
        g = math.exp(-0.25) * dw / 2
        tmpR = 0
        while (tmpR < Rm):
            n += 1
            Nc = (n - 1) % Nl
            nLayer = math.floor((n - 1) / Nl)
            nx = Nc * k
            ny = r0 + k * nLayer
            Lns = Mut(ny, ny, g, 0)
            # self inductance of current turn
            lw = lw + 2 * M_PI * ny
            M = 0
            if (n > 1):
                for j in range(n, 1, -1):
                    Jc = (j - 2) % Nl
                    jx = Jc * k
                    jLayer = math.floor((j - 2) / Nl)
                    jy = r0 + k * jLayer
                    M = M + 2 * Mut(ny, jy, nx - jx, g)

            Ltotal += Lns + M
            Resistivity = mtrl[Material["Cu"]][Rho] * 1e2
            tmpR = (Resistivity * lw * 4) / (M_PI * dw * dw)

        bTmp = (nLayer + 1) * k
        if (bTmp > c):
            break
        if (nTmp < c):
            nTmp = (nLayer + 2) * k
            result.N = Ltotal

    N1 = n - Nl
    N2 = n + Nl
    result.sec = Ltotal
    result.thd = N1
    result.fourth = N2


def getMultilayerN_Foil(D: float, w: float, t: float, ins: float, I: float, result: _CoilResult):
    D = D / 10
    w = w / 10
    t = t / 10
    ins = ins / 10
    g = math.exp(-1.5) * w
    # g = 0.2235*(w + t)
    k = ins + t
    N = 0
    r0 = (D + t) / 2
    Ltotal = 0
    while (Ltotal <= I):
        N += 1
        ny = r0 + k * (N - 1)
        Lns = Mut(ny, ny, g, 0)
        M = 0
        if (N > 1):
            for j in np.arange(N, 1, -2):
                jy = r0 + k * (j - 2)
                r = ny - jy
                gmd = math.exp(((r * r) / (w * w)) * math.log(r) + 0.5 * (1 - ((r * r) / (w * w)))
                               * math.log(w * w + r * r) + (2 * r / w) * math.atan(w / r) - 1.5)
                gmr = math.sqrt(ny*jy)
                ra = (gmd + math.sqrt(gmd * gmd + 4 * gmr * gmr)) / 2
                rb = ra - gmd
                M = M + 2 * Mut(ra, rb, 0, gmd)

        Ltotal += Lns + M

    th = k * (N - 1)
    Do = (D + 2 * th) * 10
    Length_spiral = find_actual_spiral_length(N, D, k) * 10
    Resistivity_cu = mtrl[Material["Cu"]][Rho]*1e2
    Resistivity_al = mtrl[Material["Al"]][Rho]*1e2
    Rdcc = (Resistivity_cu * Length_spiral) / (w * t) / \
        10  # resistance of the copper foil
    Rdca = (Resistivity_al * Length_spiral) / (w * t) / \
        10  # resistance of the aliminum foil
    result.N = N
    result.sec = Length_spiral/1000
    result.thd = Do
    result.fourth = Rdcc
    result.five = Rdca


def getMultilayerI_Foil(D: float, w: float, t: float, ins: float, _N: int, result: _CoilResult):
    D = D / 10
    w = w / 10
    t = t / 10
    ins = ins / 10
    g = math.exp(-1.5) * w
    k = ins + t
    r0 = (D + t) / 2
    Ltotal = 0
    for N in range(1, _N + 1):
        ny = r0 + k * (N - 1)
        Lns = Mut(ny, ny, g, 0)
        M = 0
        if (N > 1):
            for j in range(N, 1, -1):
                jy = r0 + k * (j - 2)
                r = ny - jy
                gmd = math.exp(((r * r) / (w * w)) * math.log(r) + 0.5 * (1 - ((r * r) / (w * w)))
                               * math.log(w * w + r * r) + (2 * r / w) * math.atan(w / r) - 1.5)
                gmr = math.sqrt(ny*jy)
                ra = (gmd + math.sqrt(gmd * gmd + 4 * gmr * gmr)) / 2
                rb = ra - gmd
                M = M + 2 * Mut(ra, rb, 0, gmd)

        Ltotal += Lns + M

    th = k * (_N - 1)
    Do = (D + 2 * th) * 10
    Length_spiral = find_actual_spiral_length(_N, D, k) * 10
    Resistivity_cu = mtrl[Material["Cu"]][Rho]*1e2
    Resistivity_al = mtrl[Material["Al"]][Rho]*1e2
    Rdcc = (Resistivity_cu * Length_spiral) / (w * t) / \
        10  # resistance of the copper foil
    Rdca = (Resistivity_al * Length_spiral) / (w * t) / \
        10  # resistance of the aliminum foil
    result.N = Ltotal
    result.sec = Length_spiral / 1000
    result.thd = Do
    result.fourth = Rdcc
    result.five = Rdca


def getFerriteN(L: float, Do: float, Di: float, h: float, dw: float, mu: float, Ch: float, result: _CoilResult):
    w = 0
    cr = Ch / M_SQRT2  # Chamfer radius
    # correction factor for the chamfer
    k = 0.8584 * math.pow(cr, 2) / (h * (Do - Di) / 2)
    he = h * (1 - k)  # correction ΣA/l with chamfer by correcting h
    w = 100 * math.sqrt(L / (2 * he * mu * math.log(Do / Di)))
    D1t = Do
    D2t = Di
    ht = h
    Lt = 0
    w1 = 0
    while (1):
        wt = M_PI * D2t / dw
        # Number of turns that fit into internal diameter
        w1 = w1 + wt  # Increasing Number of turns counter
        if (w1 >= w):
            Lt = Lt + (w - (w1 - wt)) * (D1t - D2t + ht * 2)
            # Wire length of this layer with previos counter
            break

        Lt = Lt + wt * (D1t - D2t + ht * 2)  # Wire length of this layer
        Lt = Lt + Lt * 0.03
        D1t = D1t + dw  # Increasing coil size
        D2t = D2t - dw
        ht = ht + 2 * dw
        if (D2t <= dw):
            Lt = -100
            break
        if not (w1 < w):
            break
    al = round(0.2 * he * mu * math.log(Do / Di))
    lw = 0.001 * Lt
    result.N = w
    result.sec = lw
    result.thd = al


def getFerriteI(N: float, Do: float, Di: float, h: float, mu: float, Ch: float, result: _CoilResult) -> float:
    cr = Ch / M_SQRT2  # Chamfer radius
    # correction factor for the chamfer
    k = 0.8584 * math.pow(cr, 2) / (h * (Do - Di) / 2)
    he = h * (1 - k)  # correction ΣA/l with chamfer by correcting h
    al = round(0.2 * he * mu * math.log(Do / Di))
    result.thd = al
    return 2e-04 * mu * he * N * N * math.log(Do / Di)


def getPCB_N(I: float, D: float, d: float, ratio: float, layout: int, result: _CoilResult):

    N = 0.5
    s = 0
    W = 0
    iTmp = 0
    while (iTmp < I):
        N = N + 0.01
        s = (D - d) / (2 * N)
        W = s * ratio
        iTmp = getPCB_I(N, d, s, layout, result)

    if (s < 0):
        N = 0

    result.N = N
    result.sec = s
    result.thd = W


def getPCB_I(N: float, _d: float, _s: float, layout: int, result: _CoilResult) -> float:

    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0

    if layout == 0:
        c1 = 1.27
        c2 = 2.07
        c3 = 0.18
        c4 = 0.13
    elif layout == 1:
        c1 = 1.0
        c2 = 2.46
        c3 = 0.0
        c4 = 0.2

    d = _d * 1e3
    s = _s * 1e3
    D = d + 2 * s * N
    Davg = (D + d) / 2
    fi = (D - d) / (D + d)
    I = mu0 * N * N * Davg * c1 * 0.5 * \
        (math.log(c2 / fi) + c3 * fi + c4 * fi * fi)

    result.five = D / 1e3
    return (I)


def log_GMD2(s: float, w: float, h: float) -> float:
    return math.log(w + h) + math.log(s / (2 * w)) - (-1.46 * w / h + 1.45) / (2.14 * w / h + 1.0)


def rectPCBSpiral(N: int, a: float, b: float, s: float, w: float, h: float) -> float:
    # /*
    #     "Inductance Formula for Rectangular Planar Spiral Inductors with Rectangular Conductor Cross Section",  H. A. Aebischer 2020
    #     https://www.researchgate.net/publication/339137261
    #     input in millimeters, output in microhenry
    # */
    try:
        a /= 1000.0
        b /= 1000.0
        s /= 1000.0
        w /= 1000.0
        h /= 1000.0
        log_GMD1 = math.log(w + h) - 3.0 / 2.0
        GMD1 = 0.2235 * (w + h)
        AMD1 = GMD1
        AMSD1sq = (1.0 / 6.0) * (w * w + h * h)
        log_GMD_L = N * (log_GMD1)
        for k in range(1, N):
            log_GMD_L += 2 * (N - k) * log_GMD2(k * s, w, h)

        log_GMD_L /= N * N
        AMSD_L = N * AMSD1sq
        for k in range(1, N):
            AMSD_L += 2.0 * (N - k) * math.pow((k * s), 2)

        AMSD_L /= N * N
        AMD_L = N * AMD1
        for k in range(1, N):
            AMD_L += 2.0 * (N - k) * math.exp(log_GMD2(k * s, w, h))

        AMD_L /= N * N
        log_GMD_a = 0
        for ks in range(-N + 1, N):
            log_GMD_a += (N - abs(ks)) * math.log(a + ks * s)

        log_GMD_a /= (N * N)
        log_GMD_b = 0
        for ks in range(-N + 1, N):
            log_GMD_b += (N - abs(ks)) * math.log(b + ks * s)

        log_GMD_b /= N * N
        AMSD_a = 0
        for ks in range(-N + 1, N):
            AMSD_a += (N - abs(ks)) * math.pow((a + ks * s), 2)

        AMSD_a /= N * N
        AMSD_b = 0
        for ks in range(-N + 1, N):
            AMSD_b += (N - abs(ks)) * math.pow((b + ks * s), 2)

        AMSD_b /= N * N
        AMD_a = 0
        for ks in range(-N + 1, N):
            AMD_a += (N - abs(ks)) * (a + ks * s)

        AMD_a /= N * N
        AMD_b = 0
        for ks in range(-N + 1, N):
            AMD_b += (N - abs(ks)) * (b + ks * s)

        AMD_b /= N * N
        La = 0.2 * (a * math.log(a + math.sqrt(a * a + AMSD_L)) -
                    a * log_GMD_L - math.sqrt(a * a + AMSD_L) + AMD_L)
        Lb = 0.2 * (b * math.log(b + math.sqrt(b * b + AMSD_L)) -
                    b * log_GMD_L - math.sqrt(b * b + AMSD_L) + AMD_L)
        Ma = 0.2 * (a * math.log(a + math.sqrt(a * a + AMSD_b)) -
                    a * log_GMD_b - math.sqrt(a * a + AMSD_b) + AMD_b)
        Mb = 0.2 * (b * math.log(b + math.sqrt(b * b + AMSD_a)) -
                    b * log_GMD_a - math.sqrt(b * b + AMSD_a) + AMD_a)
        L = 2 * N * N * (La + Lb - (Ma + Mb))
        return L
    except:
        # catch (const std::exception&):
        return 0


def getPCB_RectI(N: int, A: float, B: float, s: float, w: float, th: float, result: _CoilResult) -> float:

    if (N < 2):
        return 0
    if (B > A):
        B_neu = A
        A = B
        B = B_neu

    a = A - (N - 1) * s
    b = B - (N - 1) * s
    rho = ((N - 1) * s + w) / (B - (N - 1) * s)

    if N == 2:
        if (rho > 0.36001):
            return 0

    elif N >= 3 and N <= 7:
        if (rho > 0.52001):
            return 0

    elif N >= 8 and N <= 12:
        if (rho > 0.78001):
            return 0

    elif N >= 13 and N <= 20:
        if (rho > 0.86001):
            return 0

    if (N >= 21):
        if (rho > ((N - 1.0) / (N + 1.0))):
            return 0

    result.five = (A - (N - 1) * s * 2)
    return rectPCBSpiral(N, a, b, s, w, th)


def getPCB_RectN(I: float, A: float, B: float, _a: float, th: float, ratio: float, result: _CoilResult):

    a = A - (A - _a) / 2.0
    iTmp = 0
    s = 0
    w = 0
    result.N = 0
    result.sec = 0
    result.thd = 0
    for N in range(2, 10000):
        s = (A - _a) / (N - 1) / 2
        w = s * ratio
        b = B - (N - 1) * s
        iTmp = rectPCBSpiral(N, a, b, s, w, th)
        if ((abs(iTmp - I) / I) < 0.05):
            result.N = N
            result.sec = s
            result.thd = w
            break


def getSpiralN(I: float, Di: float, dw: float, s: float, result: _CoilResult):
    Di = Di / 10
    dw = dw / 10
    s = s / 10
    g = math.exp(-0.25) * dw / 2
    k = s + dw
    N = 0
    r0 = (Di + dw) / 2
    w = 0
    Ltotal = 0
    while (Ltotal < I):
        N += 1
        ny = r0 + k * (N - 1)
        Lns = Mut(ny, ny, g, 0)
        M = 0
        if (N > 1):
            for j in range(N, 1, -1):
                jy = r0 + k * (j - 2)
                M = M + 2 * Mut(ny, jy, 0, g)

        Ltotal += Lns + M

    w = k * (N - 1)
    Do = (Di + 2 * w) * 10
    Length_spiral = find_actual_spiral_length(N, Di, k) * 10
    result.N = N
    result.sec = Length_spiral/1000
    result.thd = Do


def getSpiralI(Do: float, Di: float, dw: float, _N: int, result: _CoilResult):

    Di = Di / 10
    Do = Do / 10
    dw = dw / 10
    g = math.exp(-0.25) * dw / 2
    w = (Do - Di) / 2
    k = w / (_N - 1)
    r0 = (Di + dw) / 2
    Ltotal = 0
    for N in range(1, _N):
        ny = r0 + k * (N - 1)
        Lns = Mut(ny, ny, g, 0)
        M = 0
        if (N > 1):
            for j in range(N, 1, -1):
                jy = r0 + k * (j - 2)
                M = M + 2 * Mut(ny, jy, 0, g)

        Ltotal += Lns + M

    Length_spiral = find_actual_spiral_length(_N, Di, k) * 10
    result.N = Ltotal
    result.sec = Length_spiral / 1000


def CalcLC0(L: float, C: float) -> float:
    f = 1e3 / (2 * M_PI * math.sqrt(L * C))
    return f


def CalcLC1(C: float, f: float) -> float:
    L1 = 1e3 / (2 * M_PI * f)
    L = math.pow(L1, 2) / C
    return L


def CalcLC2(L: float, f: float) -> float:
    C1 = 1e3 / (2 * M_PI * f)
    C = math.pow(C1, 2) / L
    return C


def CalcLC3(Zo: float, f: float, result: _CoilResult):
    result.N = 1e6 / (2 * M_PI * f * Zo)
    result.sec = Zo / (2 * M_PI * f)


def findToroidPemeability(N: float, I: float, Do: float, Di: float, h: float, Ch: float, result: _CoilResult):
    cr = Ch / M_SQRT2  # Chamfer radius
    # correction factor for the chamfer
    k = 0.8584 * math.pow(cr, 2) / (h * (Do - Di) / 2)
    he = h * (1 - k)  # correction ΣA/l with chamfer by correcting h
    m = math.ceil(10000 * I / (2 * N * N * he * math.log(Do / Di)))
    al = 1000 * I / (N * N)
    result.N = m
    result.sec = al


def findFerriteRodN(I: float, Lr: float, Dr: float, mu: float, dc: float, s: float, dw: float, p: float, result: _CoilResult):
    # Based on "The Inductance of Ferrite Rod Antennas Issue" by Alan Payne
    # [10.1][10.2] http://g3rbj.co.uk/wp-content/uploads/2014/06/Web-The-Inductance-of-Ferrite-Rod-Antennas-issue-3.pdf

    x2 = 2 * s / Lr
    N = 0
    Dk = dc + dw
    iTmp = 0
    Lf_Lair = 0
    lc = 0
    # A.V. Cainov's regression analysis algorithm calculates correction factor on shifting the coil from the center of the rod
    kx = -440.9943706 * math.pow(x2, 8) + 1318.707293 * math.pow(x2, 7) - 1604.5491034 * math.pow(x2, 6) + 1021.078226 * math.pow(
        x2, 5) - 363.8218957 * math.pow(x2, 4) + 71.6178135 * math.pow(x2, 3) - 7.6027344 * math.pow(x2, 2) + 0.3013663 * x2 + 0.995
    if ((kx < 0) or (kx > 1)):
        return

    e0 = 8.8542e-12
    dLp = -1e-4 * Dk * ((p / dw) - 1) * ((12 - (p / dw)) / 4)
    while (iTmp < I):
        N += 1
        lc = N * p
        k = lc / Dk
        # The optimized version of Wheeler's Continuous Inductance formula for the one-layer coil corrected by Robert Weaver
        # [34] (http://electronbunker.ca/eb/CalcMethods3b.html)
        i = 0.0002 * M_PI * Dk * N * N * (math.log(1 + M_PI / (2 * k)) + 1 / (2.3004 + 3.437 * k + 1.763 *
                                                                              k * k - 0.47 / math.pow((0.755 + 1 / k), 1.44)))

        Ks = rosaKs(p / dw)
        Km = rosaKm(N)
        Lcor = 0.0002 * M_PI * Dk * N * (Ks + Km)
        Lc = i - Lcor
        if (N > 4):
            Lc = Lc + dLp * N

        if (mu > 1):
            Canf = 5e-4 * M_PI * e0 * (Lr - lc) / \
                (math.log(2 * (Lr + Dr) / Dr) - 1)
            phi_phimax = 1 / (1 + (math.pow((Lr - lc) / Dr, 1.4) / (5 * mu)))
            _k = ((phi_phimax * Canf / e0) + 2e-3 * Dr) / (2e-3 * dc)
            l_prim = lc + 0.45 * dc
            _x = 5.1 * (l_prim / dc) / (1 + 2.8 * (dc / l_prim))
            mufe = (mu - 1) * math.pow(Dr / dc, 2) + 1
            Lf_Lair = (1 + _x) / (1 / _k + _x / mufe)
            iTmp = Lc * Lf_Lair * kx
        else:
            Lf_Lair = 1
            iTmp = Lc

    result.N = N
    result.sec = Lf_Lair
    result.thd = lc


def findMeadrPCB_I(a: float, d: float, h: float, W: float, N: int, result: _CoilResult):
    # http://www.journal.ftn.kg.ac.rs/Vol_1-3/08-Stojanovic-Zivanov-Damnjanovic.pdf (The monomial equation [11])

    result.N = 0.00266 * math.pow(a, 0.0603) * math.pow(h, 0.4429) * \
        math.pow(N, 0.954) * math.pow(d, 0.606) * math.pow(W, -0.173)
    result.sec = 2 * N * d + 2 * a


def findMultiloop_I(N: float, Di: float, dw: float, dt: float, result: _CoilResult) -> float:
    # The author of the source code of this function is George Overton.
    # The source code is used as an open with the consent of the author.
    # The code is from the author's book "Inside the METAL DETECTOR" Appendix A
    # https://www.facebook.com/Inside-the-Metal-Detector-222330481232886/
    # https://www.geotech1.com

    c = math.sqrt(N) * dt
    a = (Di + c) / 2
    x = math.pow(c / 2 / a, 2)
    s1 = 0.0004 * M_PI * a
    s2 = math.pow(N, 2)
    s3 = s1 * s2
    s4 = 0.5 + x/12
    s5 = math.log(8/x)
    s6 = (s4 * s5) - 0.85 + (0.2 * x)
    ind = s3 * s6  # Inductance (microH)
    result.N = 2 * a  # Mean coil diameter (mm)
    result.sec = c  # coil thickness (mm)
    result.thd = 2e-3 * M_PI * a * N  # length of the wire (m)
    Resistivity_cu = mtrl[Material["Cu"]][Rho]*1e2
    result.fourth = (Resistivity_cu * result.thd * 100 * 4) / \
        (M_PI * dw * dw * 0.01)  # Resistance to DC (Ohm)
    return ind


def findMultiloop_N(I: float, Di: float, dw: float, dt: float, result: _CoilResult) -> float:
    tmpI = 0
    N = 1
    i = 0
    while (tmpI <= I):
        i += 1
        tmpI = findMultiloop_I(N, Di, dw, dt, result)
        N += 0.1
        if ((N > 1e7) or ((i == 1) and (tmpI > I))):
            return -1

    return N


def findRoundLoop_I(D: float, dw: float) -> float:
    #

    return 0.0002 * M_PI * D * (math.log(8 * D / dw) - 2)


def findRoundLoop_D(Ind: float, dw: float) -> float:
    tmpI = 0
    D = 2 * dw
    i = 0
    while (tmpI <= Ind):
        i += 1
        tmpI = findRoundLoop_I(D, dw)
        D += 0.01
        if ((D > 2e4) or ((i == 1) and (tmpI > Ind))):
            return -1

    return D


def findIsoIsoscelesTriangleLoop_I(_a: float, _b: float, dw: float) -> float:
    #

    c = _a / 1000
    b = _b / 1000
    r = dw / 2000

    a1 = 2 * c * math.log(2 * c / r)
    a2 = b * math.log(2 * c / r)
    a3 = 2 * (b + c) * math.asinh(b * b /
                                  (math.sqrt(4 * b * b * c * c - math.pow(b, 4))))
    a4 = 2 * c * math.asinh((2 * c * c - b * b) /
                            (math.sqrt(4 * b * b * c * c - math.pow(b, 4))))
    a5 = 2 * c + b
    return 0.2 * (a1 + a2 - a3 - a4 - a5)


def findIsoIsoscelesTriangleLoop_a(Ind: float, dw: float) -> float:
    tmpI = 0
    a = 2 * dw
    i = 0
    while (tmpI <= Ind):
        i += 1
        tmpI = findIsoIsoscelesTriangleLoop_I(a, a, dw)
        a += 0.01
        if ((a > 2e4) or ((i == 1) and (tmpI > Ind))):
            return -1

    return a


def findRectangleLoop_I(_a: float, _b: float, dw: float) -> float:
    #

    a = _a / 1000
    b = _b / 1000
    r = dw / 2000
    d = math.sqrt(a * a + b * b)

    a1 = -2 * (a + b)
    a2 = 2 * d
    a3 = b * math.log((b + d) / a)
    a4 = a * math.log((a + d) / b)
    a5 = b * math.log(2 * b / r)
    a6 = a * math.log(2 * a / r)
    return 0.4 * (a1 + a2 - a3 - a4 + a5 + a6)


def findRectangleLoop_a(Ind: float, dw: float) -> float:
    tmpI = 0
    a = 2 * dw
    i = 0
    while (tmpI <= Ind):
        i += 1
        tmpI = findRectangleLoop_I(a, a, dw)
        a += 0.01
        if ((a > 2e4) or ((i == 1) and (tmpI > Ind))):
            return -1

    return a


def findSheildedInductance(I: float, D: float, Ds: float, l: float, Hs: float) -> float:
    result = I * (1 - (D / Ds) * (D / Ds) * (D / Ds)) * \
        (1 - (l / (2 * Hs)) * (l / (2 * Hs)))
    return result


def findAirCoreRoundToroid_I(N: float, D1: float, D2: float, dw: float) -> float:
    R = 0.1 * (D1 + D2) / 4
    a = 0.05 * (((D1 - D2) / 2) + dw)
    ind = 0.01257 * N * N * (R - math.sqrt(R * R - a * a))
    return ind


def findAirCoreRoundToroid_N(Ind: float, D1: float, D2: float, dw: float) -> float:
    R = 0.1 * (D1 + D2) / 4
    a = 0.05 * (((D1 - D2) / 2) + dw)
    N = math.sqrt(Ind / (0.01257 * (R - math.sqrt(R * R - a * a))))
    return N


def findPotCore_I(N: float, D1: float, D2: float, D3: float, D4: float, h1: float, h2: float, g: float, b: float, mu: float, result: _CoilResult) -> float:
    r1 = 0.5 * D4
    r2 = 0.5 * D3
    r3 = 0.5 * D2
    r4 = 0.5 * D1
    h = 0.5 * (h1 - h2)
    l1 = h2
    l3 = h2
    k1 = 2 * b * (r4 - r3)
    A1 = M_PI * (r4 - r3) * (r4 + r3) - k1
    A3 = M_PI * (r2 - r1) * (r2 + r1)
    sum11 = l1 / A1
    sum21 = l1 / (A1 * A1)
    k2 = 1 / (1 - (2 * b / (2 * M_PI * r3)))
    sum12 = k2 * math.log(r3 / r2) / (M_PI * h)
    sum22 = k2 * (r3 - r2) / (2 * r3 * r2 * M_PI * M_PI * h * h)
    sum13 = l3 / A3
    sum23 = l3 / (A3 * A3)
    s1 = r2 - math.sqrt(0.5 * (r2 * r2 + r1 * r1))
    s2 = math.sqrt(0.5 * (r3 * r3 + r4 * r4)) - r3
    l4 = 0.25 * M_PI * (2 * s2 + h)
    l5 = 0.25 * M_PI * (2 * s1 + h)
    k4 = 1 - (2 * b / (M_PI * (r3 + r4)))
    A4 = 0.5 * k4 * M_PI * (r4 * r4 - r3 * r3 + 2 * r3 * h)
    A5 = 0.5 * M_PI * (r2 * r2 - r1 * r1 + 2 * r2 * h)
    sum14 = l4 / A4
    sum24 = l4 / (A4 * A4)
    sum15 = l5 / A5
    sum25 = l5 / (A5 * A5)
    C1 = sum11 + sum12 + sum13 + sum14 + sum15
    C2 = sum21 + sum22 + sum23 + sum24 + sum25
    le = C1 * C1 / C2
    Ae = C1 / C2
    mu_e = mu / (1 + g * mu / le)
    ind = 1000 * N * N * mu0 * mu_e / C1
    result.N = le
    result.sec = Ae
    result.thd = mu_e
    return ind


def findPotCore_N(Ind: float, D1: float, D2: float, D3: float, D4: float, h1: float, h2: float, g: float, b: float, mu: float, result: _CoilResult):
    tmpI = 0
    N = 0
    while (tmpI <= Ind):
        N += 1
        tmpI = findPotCore_I(N, D1, D2, D3, D4, h1, h2, g, b, mu, result)

    return N


def findECore_I(N: float, A: float, B: float, C: float, D: float, E: float, F: float, g: float, b: float, mu: float,
                isEI: bool, isRound: bool, result: _CoilResult) -> float:

    k = 1
    if (isRound):
        k = 1.1918
    if (isEI):
        l1 = D
        l3 = D
    else:
        l1 = 2 * D
        l3 = 2 * D

    h = B - D
    s = 0.5 * F
    p = (A - E) / 2
    q = C
    l2 = E - F
    l4 = M_PI * (p + h) / 4
    l5 = M_PI * (k * s + h) / 4
    A1 = 2 * q * p - 0.5 * M_PI * b * b
    A2 = 2 * q * h
    if (isRound):
        A3 = M_PI * s * s
    else:
        A3 = 2 * s * q
    A4 = 0.5 * (A1 + A2)
    A5 = 0.5 * (A2 + A3)
    c = getFerriteCoreMagConst(l1, l2, l3, l4, l5, A1, A2, A3, A4, A5)
    le = c.C1 * c.C1 / c.C2
    Ae = c.C1 / c.C2
    mu_e = mu / (1 + g * mu / le)
    ind = 1000 * N * N * mu0 * mu_e / c.C1
    result.N = le
    result.sec = Ae
    result.thd = mu_e
    return ind


def findECore_N(Ind: float, A: float, B: float, C: float, D: float, E: float, F: float, g: float, b: float, mu: float,
                isEI: bool, isRound: bool, result: _CoilResult):
    tmpI = 0
    N = 0
    while (tmpI <= Ind):
        N += 1
        tmpI = findECore_I(N, A, B, C, D, E, F, g, b,
                           mu, isEI, isRound, result)

    return N


def findUCore_I(N: float, A: float, B: float, C: float, D: float, E: float, F: float, s: float, mu: float, result: _CoilResult) -> float:
    l1 = 2 * D
    l3 = 2 * D
    h = B - D
    q = C
    l2 = 2 * E
    p = 0
    y = 0
    A1 = 0
    A3 = 0
    if (F <= 0):
        p = (A - E) / 2
        y = p
    else:
        p = A - E - F
        y = F

    l4 = M_PI * (p + h) / 4
    l5 = M_PI * (y + h) / 4
    if (F == 0):
        A1 = q * p
    elif (F < 0):
        A1 = 0.25 * M_PI * q * q - 0.25 * M_PI * s * s
    elif (F > 0):
        A1 = 0.25 * M_PI * q * q
    A2 = q * h
    if (F == 0):
        A3 = q * p
    elif (F < 0):
        A3 = 0.25 * M_PI * q * q - 0.25 * M_PI * s * s
    elif (F > 0):
        A3 = q * y - 0.25 * M_PI * s * s
    A4 = 0.5 * (A1 + A2)
    A5 = 0.5 * (A2 + A3)
    c = getFerriteCoreMagConst(l1, l2, l3, l4, l5, A1, A2, A3, A4, A5)
    le = c.C1 * c.C1 / c.C2
    Ae = c.C1 / c.C2
    ind = 1000 * N * N * mu0 * mu / c.C1
    result.N = le
    result.sec = Ae
    return ind


def findUCore_N(Ind: float, A: float, B: float, C: float, D: float, E: float, F: float, s: float, mu: float, result: _CoilResult):
    tmpI = 0
    N = 0
    while (tmpI <= Ind):
        N += 1
        tmpI = findUCore_I(N, A, B, C, D, E, F, s, mu, result)

    return N


def findRMCore_I(N: float, a: float, b: float, c: float, e: float, d2: float, d3: float, d4: float, h1: float,  h2, g: float, mu: float, type: int, result: _CoilResult) -> float:

    h = 0.5 * (h1 - h2)
    alpha = M_PI / 2.0
    betta = alpha - math.asin(e / d2)
    phi = M_PI / 2.0
    if (type == 2):
        phi = 2 * math.acos(e / d2)
    p = math.sqrt(2.0) * a - b
    lmin = 0.5 * (d2 - d3)
    lmax = math.sqrt(0.25 * (d2 * d2 + d3 * d3) - 0.5 *
                     d2 * d3 * math.cos(alpha - betta))
    if (type == 2):
        lmax = math.sqrt(0.25 * (d2 * d2 + d3 * d3) - 0.5 * d2 * d3 *
                         math.cos(alpha - betta)) - 0.5 * b / math.sin(0.5 * phi)
    if (type == 3):
        lmax = 0.5 * e + 0.5 * (1 - math.sin(M_PI / 4.0)*(d2 - c))
    f = (lmax + lmin) / (2.0 * lmin)
    A7 = 0.25 * (0.5 * betta * d2 * d2 + 0.5 * e * e * math.tan(betta) -
                 0.5 * e * e * math.tan(alpha - 0.5 * phi)-0.25 * M_PI * d3 * d3)
    if (type == 1):
        A7 = 0.25 * (0.5 * betta * d2 * d2 + 0.5 * d2 * d3 * math.sin(alpha - betta) +
                     0.5 * math.pow(c - d3, 2) * math.tan(0.5 * phi) - 0.25 * M_PI * d3 * d3)
    if (type == 2):
        A7 = 0.25 * (0.5 * betta * d2 * d2 - 0.25 * M_PI * d3 * d3 + 0.5 * (b *
                     b - e * e) * math.tan(alpha - 0.5 * phi) + 0.5 * e * e * math.tan(betta))
    if (type == 3):
        A7 = 0.25 * (0.5 * betta * d2 * d2 - 0.5 * phi * d3 *
                     d3 + 0.5 * c * c * math.tan(alpha - betta))
    A8 = (M_PI / 16) * (d2 * d2 - d3 * d3)
    D = A7 / A8

    l1 = h2
    A1 = 0.5 * a * a * (1 + math.tan(betta - M_PI / 4.0)) - \
        0.5 * betta * d2 * d2 - 0.5 * p * p
    sum11 = l1 / A1
    sum12 = l1 / (A1 * A1)

    sum21 = (math.log(d2 / d3) * f) / (M_PI * D * h)
    sum22 = (f * (1.0 / d3 - 1.0 / d2)) / (math.pow(M_PI * D * h, 2))

    l3 = h2
    A3 = 0.25*M_PI * (d3 * d3 - d4 * d4)
    sum31 = l3 / A3
    sum32 = l3 / (A3 * A3)

    l4 = 0.25 * M_PI * (h + 0.5 * a - 0.5 * d2)
    A4 = 0.5 * (A1 + 2.0 * betta * d2 * h)
    sum41 = l4 / A4
    sum42 = l4 / (A4 * A4)

    l5 = 0.25 * M_PI * (d3 + h - math.sqrt(0.5 * (d3 * d3 + d4 * d4)))
    A5 = 0.5 * (0.25 * M_PI * (d3 * d3 - d4 * d4) + 2.0 * alpha * d3 * h)
    sum51 = l5 / A5
    sum52 = l5 / (A5 * A5)

    C1 = sum11 + sum21 + sum31 + sum41 + sum51
    C2 = sum12 + sum22 + sum32 + sum42 + sum52
    le = C1 * C1 / C2
    Ae = C1 / C2
    mu_e = mu / (1 + g * mu / le)
    ind = 1000 * N * N * mu0 * mu_e / C1
    result.N = le
    result.sec = Ae
    result.thd = mu_e
    return ind


def findRMCore_N(Ind: float, a: float, b: float, c: float, e: float, d2: float, d3: float, d4: float, h1: float,  h2, g: float, mu: float, type: int, result: _CoilResult):
    tmpI = 0
    N = 0
    while (tmpI <= Ind):
        N += 1
        tmpI = findRMCore_I(N, a, b, c, e, d2, d3, d4,
                            h1, h2, g, mu, type, result)

    return N


def findBrooksCoil(I: float, d: float, pa: float, pr: float):
    a = 0.0025491
    _I = 0
    ha = pa * d
    hr = pr * d
    N = 0

    N += 0.01
    discr = math.pow(ha + d, 2) - 4 * (ha * d - hr * ha * N)
    c = 0.5 * (ha + d + math.sqrt(discr))
    _I = a * c * N * N
    while (I > _I):
        N += 0.01
        discr = math.pow(ha + d, 2) - 4 * (ha * d - hr * ha * N)
        c = 0.5 * (ha + d + math.sqrt(discr))
        _I = a * c * N * N

    Nc = (c / ha - 1)   # number of turns in the layer
    nLayer = math.ceil(N / Nc)     # number of layers
    # Calculation of the wire length in all layers except the last
    if (nLayer > 1):
        for j in range(nLayer - 2 + 1):
            lengthWire += math.sqrt(ha * ha +
                                    math.pow(M_PI * (2 * c + d + 2 * hr * j), 2))

    else:
        lengthWire = 0
    lengthWire = Nc * lengthWire
    # Calculating the wire length in the last layer
    lengthWireLastLayer = (math.sqrt(
        ha * ha + math.pow(M_PI * (2 * c + d + 2 * hr * (nLayer - 1)), 2))) * (N - Nc * (nLayer - 1))
    lengthWire += lengthWireLastLayer
    lengthWire = lengthWire / 1000
    Resistivity = mtrl[Material["Cu"]][Rho]*1e6
    # Resistance of the wire to DC [Ohm]
    DCR = (Resistivity * lengthWire * 4) / (M_PI * d * d)
    massWire = 2.225 * M_PI * d * d * \
        lengthWire          # Weight of the wire [g]

    return N, nLayer, Nc, c, lengthWire, massWire, DCR


def findPadderCapacitance(Ct: float, Cv_low: float, Cv_high: float, Cstray: float, cap_ratio: float) -> float:
    chp = Cv_high + Ct
    clp = Cv_low + Ct
    beta = (cap_ratio - 1) * Cstray
    a = cap_ratio * clp - chp + beta
    b = (cap_ratio - 1) * clp * chp + beta * (clp + chp)
    c = beta * clp * chp
    return (-b - math.sqrt(b * b - 4 * a * c)) / (2.0 * a)


def findTrimmerCapacitance(Cp: float, Cv_low: float, Cv_high: float, Cstray: float, cap_ratio: float) -> float:
    Ct = 0.0
    if (Cp < 0.00001):
        Ct = (Cv_high + Cstray - cap_ratio *
              (Cv_low + Cstray)) / (cap_ratio - 1)
    else:
        k0 = (1 - cap_ratio) * Cstray
        k1 = k0 * Cp * Cp
        k2 = k0 * Cp + Cp * Cp
        k3 = k0 * Cp - cap_ratio * Cp * Cp
        k4 = k0 + (1 - cap_ratio) * Cp
        a = k4
        b = k2 + k3 + k4 * (Cv_low + Cv_high)
        c = k1 + k2 * Cv_high + k3 * Cv_low + k4 * Cv_low * Cv_high
        Ct = (-b - math.sqrt(b * b - 4 * a * c)) / (2.0 * a)

    return Ct
