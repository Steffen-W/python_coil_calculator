

# from matplotlib import patches
# import matplotlib.pyplot as plt

# the library must be compiled
# g++ -shared -o coil64_lib.so -fPIC -I /usr/include/python3.10 -lpython3.10 *.cpp
import coil64_lib

from enum import Enum
import numpy as np


class Material(Enum):
    Al = 0
    Cu = 1
    Ag = 2
    Ti = 3


class layoutType(Enum):
    Square = 0
    Spiral = 1
    Rectangular = 2


def calc_Onelayer_cw(I=10.0, f=1.0, D=15.0, l=10.0, mt=Material.Cu):
    """One layer close-winding coil
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil1.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        f (float): Frequency f (MHz). Defaults to 1.0.
        D (float): Former diameter (mm). Defaults to 15.0.
        l (float): Winding length (mm). Defaults to 10.0.
        mt (_type_, optional): Material. Defaults to Material.Cu.

    Returns:
        _type_: result
    """

    result = coil64_lib.getOneLayerN_byWindingLength(D, l, I)
    od = coil64_lib.odCalc(result["dw"])
    Dk = D + od

    Cs = coil64_lib.find_Cs(od, Dk, od * result["N"])
    result2 = coil64_lib.solve_Qr(
        I, Dk, od, result["dw"], f, result["N"], Cs, mt.value)
    srf = coil64_lib.findSRF(od * result["N"], Dk, result["lw"])

    print("Number of turns of the coil {:.3f}".format(result["N"]))
    print("Length of wire without leads lw = {:.3f} mm".format(
        1000*result["lw"]))
    print("Self capacitance Cs = {:.3f} pF".format(Cs))
    print("Coil self-resonance frequency Fsr = {:.3f} MHz".format(srf))
    print("Coil constructive Q-factor = {:.3f}".format(result2["R_ind / Rac"]))
    print("Loss resistance ESR = {:.3f} Ohm".format(result2["Rac"]))
    result["R_ind / Rac"] = result2["R_ind / Rac"]
    result["Rac"] = result2["Rac"]
    result["Cs"] = Cs
    result["srf"] = srf
    return result


# def calc_Onelayer_q():  # One layer coil on a polygonal former
#     # https://coil32.net/one-layer-air-core-coil.html
#     # /images/res/Coil3.png

#     I = 10  # Inductance L
#     f = 1  # Frequency f
#     D = 10  # Former diameter
#     d = 0.1  # Wire diameter
#     k = 0.11  # Wire diameter with insulation
#     p = 1.5  # Winding pitch
#     n = 6  # Number of sides of the former

#     mt = resolve_q.Material.Cu
#     accuracy = 1

#     # number of turns
#     result.N = resolves.getOneLayerN_Poligonal(I, D, d, p, n, result, accuracy)
#     result.fourth = resolve_srf_cs.find_Cs(
#         p, result.seven, p * result.N)    # self-capacitance
#     # self-resonance frequency    Q-factor
#     result.five = resolve_srf_cs.findSRF(
#         p * result.N, result.seven, result.thd)
#     result.six = resolve_q.solve_Qr(
#         I, result.seven, p, d, f, result.N, result.fourth, mt, result)

#     print("Number of turns of the coil N = {:.3f}".format(result.N))# from matplotlib import patches
# import matplotlib.pyplot as plt
#     print("Length of wire without leads lw = {:.3f} mm".format(1000*result.thd))
#     print("Length of winding l = {:.3f} mm".format(result.sec))
#     # print("DC resistance of the coil Rdc = {:.3f} Ohm".format(result.))
#     # print("Reactance of the coil X = {:.3f} Ohm".format(result.))
#     print("Self capacitance Cs = {:.3f} pF".format(result.fourth))
#     print("Coil self-resonance frequency Fsr = {:.3f} MHz".format(result.five))
#     print("Coil constructive Q-factor Q = {:.3f}".format(result.six))
#     print("Loss resistance ESR = {:.3f} Ohm".format(result.seven))

def calc_Onelayer_p(I=10, f=1.0, D=10.0, w=1.0, t=0.1, i=0.1, p=1.5, mt=Material.Cu):
    """One layer coil with rect wire
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil2_square.png

    Args:
        I (int): Inductance L (uH). Defaults to 10.
        f (float): Frequency f (MHz). Defaults to 1.0.
        D (float): Former diameter (mm). Defaults to 10.0.
        w (float): Wire width (mm). Defaults to 1.0.
        t (float): Wire thickness (mm). Defaults to 0.1.
        i (float, optional): Insulation thickness (mm). Defaults to 0.1.
        p (float, optional): Winding pitch (mm). Defaults to 1.5.
        mt (_type_, optional): _description_. Defaults to Material.Cu.

    Returns:
        _type_: result
    """

    Dk = D
    t = i

    result = coil64_lib.getOneLayerN_withRectWire(
        D, w, t, p, I)
    result["Cs"] = coil64_lib.find_Cs(p, D, p * result["N"])
    result2 = coil64_lib.solve_Qc(
        I, D, p, w, t, f, result["N"], result["Cs"], mt.value)
    result["SRF"] = coil64_lib.findSRF(p * result["N"], D, result["lw"])

    result["Q"] = result2["Q"]
    result["Rac"] = result2["Rac"]

    print("Number of turns of the coil N = {:.3f}".format(result["N"]))
    print("Length of wire without leads lw = {:.3f} mm".format(
        1000*result["lw"]))
    # print("Length of winding l = 115.43 mm".format(result.))
    # print("DC resistance of the coil Rdc = 0.418 Ohm".format(result.))
    # print("Reactance of the coil X = 31.416 Ohm".format(result.))

    print("Self capacitance Cs = 1.418 pF".format(result["Cs"]))
    print(
        "Coil self-resonance frequency Fsr = {:.3f} MHz".format(result["SRF"]))
    print("Coil constructive Q-factor Q = {:.3f}".format(result["Q"]))
    print("Loss resistance ESR = {:.3f} Ohm".format(result["Rac"]))
    return result


def calc_L_Onelayer(N=15, f=1.0, D=15.0, d=0.1, k=0.11, p=1.5, mt=Material.Cu):
    """calc_L_Onelayer
    One layer coil with round wire
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil2.png

    Args:
        N (int): Number of turns. Defaults to 15.
        f (float): Frequency. Defaults to 1.0.
        D (float): Former diameter. Defaults to 15.0.
        d (float): Wire diameter. Defaults to 0.1.
        k (float, optional): Wire diameter with insulation. Defaults to 0.11. # TODO to check why unused
        p (float): Winding pitch. Defaults to 1.5.
        mt (_type_, optional): Material. Defaults to Material.Cu.

    Returns:
        _type_: result
    """

    result = coil64_lib.getOneLayerI_withRoundWire(D, d, p, N)
    result["Cs"] = coil64_lib.find_Cs(p, D, p * N)  # self-capacitance
    result["SRF"] = coil64_lib.findSRF(
        p * N, D, result["lw"])  # self-resonance frequency
    result2 = coil64_lib.solve_Qr(
        result["L"], D, p, d, f, N, result["Cs"], mt.value)  # Q-factor
    result["Rac"] = result2["Rac"]
    result["Q"] = result2["R_ind / Rac"]

    print("Inductance L = {:.3f} uH".format(result["L"]))
    print("Length of wire without leads lw = {:.3f} mm".format(
        1000*result["lw"]))
    print("Self capacitance Cs = {:.3f} pF".format(result["Cs"]))
    print(
        "Coil self-resonance frequency Fsr = {:.3f} MHz".format(result["SRF"]))
    print("Coil constructive Q-factor Q = {:.3f}".format(result["Q"]))
    print("Loss resistance ESR = {:.3f} Ohm".format(result["Rac"]))
    return result


def calc_Multilayer(I=20.0, D=8.0, lk=5.0, d=0.1, k=0.109):
    """Multilayer coil
    https://coil32.net/multi-layer-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil4.png

    Args:
        I (float): Inductance L (uH). Defaults to 20.0.
        D (float): Former diameter (mm). Defaults to 8.0.
        lk (float): Winding length (mm). Defaults to 5.0.
        d (float): Wire diameter (mm). Defaults to 0.1.
        k (float): Wire diameter with insulation (mm). Defaults to 0.109.

    Returns:
        _type_: result
    """

    result = coil64_lib.getMultiLayerN(I, D, d, k, lk, 0, -1)

    print("Number of turns of the coil N = {:.3f}".format(
        result['Number turns']))  # TODO: can be wrong!
    print("Thickness of the coil c = {:.3f} mm".format(result['Thickness']))
    print("Length of wire without leads lw = {:.3f} mm".format(
        1000*result['Length']))
    print("DC resistance of the coil Rdc = {:.3f} Ohm".format(result['R_DC']))
    print("Number of layers Nl = {:.3f}".format(result['Number layers']))
    return result


def calc_Multilayer_p(I=10.0, D=15.0, lk=8.0, d=0.1, k=0.11, g=0.1, Ng=5.0):
    """Multilayer coil with insulated pads
    https://coil32.net/multi-layer-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil4-0.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        D (float): Former diameter (mm). Defaults to 15.0.
        lk (float): Winding length (mm). Defaults to 8.0.
        d (float): Wire diameter (mm). Defaults to 0.1.
        k (float): Wire diameter with insulation (mm). Defaults to 0.11.
        g (float, optional): Insulation thickness (mm). Defaults to 0.1.
        Ng (float, optional): Layers number beetween insulating pads (-). Defaults to 5.0.

    Returns:
        _type_: result
    """

    result = coil64_lib.getMultiLayerN(I, D, d, k, lk, g, Ng)

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Thickness of the coil c = {:.3f} mm".format(result["Thickness"]))
    print("Length of wire without leads lw = {:.3f} mm".format(
        1000*result["Length"]))
    print("DC resistance of the coil Rdc = {:.3f} Ohm".format(result["R_DC"]))
    print("Number of layers Nl = {:.3f}".format(result["Number layers"]))
    print("Number of interlayers Ng = {:.3f}".format(result["Ng"]))
    return result


def calc_Multilayer_r(I=10.0, a=5.0, b=2.0, l=8.0, d=0.1, k=0.11):
    """Multilayer coil on a rectangular former
    https://coil32.net/multilayer-rectangular.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil4_square.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        a (float): Former width (mm). Defaults to 5.0.
        b (float): Former height (mm). Defaults to 2.0.
        l (float): Winding length (mm). Defaults to 8.0.
        d (float): Wire diameter (mm). Defaults to 0.1.
        k (float): Wire diameter with insulation (mm). Defaults to 0.11.

    Returns:
        _type_: result
    """

    result = coil64_lib.getMultiLayerN_rectFormer(I, a, b, l, d, k)

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Number of layers Nl = {:.3f}".format(result["Number layers"]))
    print("Thickness of the coil c = {:.3f} mm".format(result["thickness"]))
    print("Length of wire without leads lw = {:.3f} mm".format(1000 *
                                                               result["Length wire"]))
    print("DC resistance of the coil Rdc = {:.3f} Ohm".format(result["Rdc"]))
    return result


def calc_Multilayer_f(I=10.0, D=15.0, w=12.0, t=0.1, g=0.01):
    """Multilayer foil-wound coil
    https://coil32.net/foil-wound-coil-calculation.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil11.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        D (float): Former diameter (mm). Defaults to 15.0.
        w (float): Foil width (mm). Defaults to 12.0.
        t (float): Foil thickness (mm). Defaults to 0.1.
        g (float): Insulation thickness (mm). Defaults to 0.01.

    Returns:
        _type_: result
    """

    ins = g
    result = coil64_lib.getMultilayerN_Foil(D, w, t, ins, I)

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Outside diameter Do = {:.3f} mm".format(result["Do"]))
    print("Length of the foil lf = {:.3f} mm".format(
        1000*result["Length spiral"]))
    print("DC resistance of the coil Rdc = {:.3f} Ohm (Copper)".format(
        result["Rdcc"]))
    print("DC resistance of the coil Rdc = {:.3f} Ohm (Aluminum)".format(
        result["Rdca"]))
    return result


def calc_L_Multilayer_f(N=10, D=15.0, w=10.0, t=0.1, g=0.01):
    """Multilayer foil-wound coil
    https://coil32.net/foil-wound-coil-calculation.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil11.png

    Args:
        N (int): Number of turns. Defaults to 10.
        D (float): Former diameter (mm). Defaults to 15.0.
        w (float): Foil width (mm). Defaults to 10.0.
        t (float, optional): Foil thickness (mm). Defaults to 0.1.
        g (float, optional): Insulation thickness (mm). Defaults to 0.01.

    Returns:
        _type_: result
    """

    ins = g
    result = coil64_lib.getMultilayerI_Foil(D, w, t, ins, N)

    print("Inductance L = {:.3f} uH".format(result["L"]))
    print("Length of the foil lf = {:.3f} mm".format(1000*result["Length"]))
    print("Outside diameter Do = {:.3f} mm".format(result["Do"]))
    print("DC resistance of the coil Rdc = {:.3f} Ohm (Copper)".format(
        result["R_DC"]))
    print("DC resistance of the coil Rdc = {:.3f} Ohm (Aluminum)".format(
        result["R_AC"]))
    return result


def calc_FerrToroid(I=10.0, OD=30.0, ID=10.0, h=5.0, d=0.0, mu=10.0, C=0.0):
    """calc_FerrToroid
    https://coil32.net/ferrite-toroid-core.html
    https://github.com/radioacoustick/Coil64/tree/master/res/T-core.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        OD (float): Outside diameter (mm). Defaults to 30.0.
        ID (float): Inside diameter (mm). Defaults to 10.0.
        h (float): Core height (mm). Defaults to 5.0.
        d (float): Wire diameter (mm). Defaults to 0.0.
        mu (float): Init magnetic permeability (-). Defaults to 10.0.
        C (float): Chamfer (mm). Defaults to 0.0.

    Returns:
        _type_: result
    """

    result = coil64_lib.getFerriteN(I, OD, ID, h, d, mu, C)

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Length of wire without leads lw = {:.3f} mm".format(1000 *
                                                               result["Length wire"]))
    print("AL = {:.3f} nH/N2".format(result["Al"]))
    return result


def calc_L_FerrToroid(N=5, OD=20.0, ID=10.0, h=5.0, mu=10.0, C=0.0):
    """calc_L_FerrToroid
    https://coil32.net/ferrite-toroid-core.html
    https://github.com/radioacoustick/Coil64/tree/master/res/T-core.png

    Args:
        N (int, optional): Number of turns. Defaults to 5.
        OD (float, optional): Outside diameter (mm). Defaults to 20.0.
        ID (float, optional): Inside diameter (mm). Defaults to 10.0.
        h (float, optional): Core height (mm). Defaults to 5.0.
        mu (float, optional): Init magnetic permeability (mm). Defaults to 10.0.
        C (float, optional): Chamferd (mm). Defaults to 0.0.

    Returns:
        _type_: result
    """

    result = coil64_lib.getFerriteI(N, OD, ID, h, mu, C)

    print("Inductance L = {:.3f} uH".format(result["L"]))
    print("A_L = {:.3f} nH/N^2".format(result["Al"]))
    return result


def calc_PCB_coil_Square(I=10.0, f=1.0, D=20.0, d=5.0, t=0.1, ratio=0.6):
    """calc_PCB_coil_Square
    https://coil32.net/pcb-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil8.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        f (float): Frequency f (MHz). Defaults to 1.0.
        D (float): Outside diameter (mm). Defaults to 20.0.
        d (float): Inside diameter (mm). Defaults to 5.0.
        t (float, optional): PCB trace thickness (mm). Defaults to 0.1.
        ratio (float, optional): wire width / wire pitch (mm). Defaults to 0.6.

    Returns:
        _type_: result
    """

    # 0: Square; 1: Spiral; 2: Rectangular
    layoutPCB = layoutType.Square

    result = coil64_lib.calc_getPCB_N(I, D, d, ratio, layoutPCB.value)

    if ((result["Winding pitch"] != 0) and (result["Width"] != 0)):
        result["Q"] = coil64_lib.solve_Qpcb(int(result["Number turns"]), I, D, d,
                                            result["Width"], t, result["Winding pitch"],
                                            f, layoutPCB.value)
    else:
        result["Q"] = 0

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Winding pitch s = {:.3f} mm".format(result["Winding pitch"]))
    print("Width of a PCB trace W = {:.3f} mm".format(result["Width"]))
    print("Coil constructive Q-factor Q ≈ {:.3f}".format(result["Q"]))
    return result


def calc_L_PCB_coil_Square(N=5, f=1.0, d=10.0, s=1.0, W=0.5, t=0.1):
    """calc_L_PCB_coil
    https://coil32.net/pcb-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil8.png

    Args:
        N (int): Number of turns. Defaults to 5
        f (float): Frequency f (MHz). Defaults to 1.0.
        d (float): Inside diameter (mm). Defaults to 10.0.
        s (float): Winding pitch (mm). Defaults to 1.0.
        W (float): Width of a PCB trace (mm). Defaults to 0.5.
        t (float, optional): PCB trace thickness (mm). Defaults to 0.1.

    Returns:
        _type_: result
    """

    # 0: Square; 1: Square; 2: Rectangular
    layoutPCB = layoutType.Square

    result = coil64_lib.getPCB_I(N, d, s, layoutPCB.value)
    A = d + 2 * s * N
    B = d
    result["Q"] = coil64_lib.solve_Qpcb(
        N, result["L"], A, B, W, t, s, f, layoutPCB.value)

    print("Inductance L (uH) = {:.3f} uH".format(result["L"]))
    print("Outside diameter D = {:.3f} mm".format(result["Do"]))
    print("Coil constructive Q-factor Q ≈ {:.3f}".format(result["Q"]))
    return result


def calc_PCB_coil_Spiral(I=10.0, f=1.0, D=20.0, d=5.0, t=0.1, ratio=0.6):
    """calc_PCB_coil_Spiral
    https://coil32.net/pcb-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil9.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        f (float): Frequency f (MHz). Defaults to 1.0.
        D (float): Outside diameter (mm). Defaults to 20.0.
        d (float): Inside diameter (mm). Defaults to 5.0.
        t (float, optional): PCB trace thickness (mm). Defaults to 0.1.
        ratio (float, optional): wire width / wire pitch (mm). Defaults to 0.6.

    Returns:
        _type_: result
    """

    layoutPCB = layoutType.Spiral

    result = coil64_lib.calc_getPCB_N(I, D, d, ratio, layoutPCB.value)

    if ((result["Winding pitch"] != 0) and (result["Width"] != 0)):
        result["Q"] = coil64_lib.solve_Qpcb(int(result["Number turns"]), I, D, d,
                                            result["Width"], t, result["Winding pitch"],
                                            f, layoutPCB.value)
    else:
        result["Q"] = 0

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Winding pitch s = {:.3f} mm".format(result["Winding pitch"]))
    print("Width of a PCB trace W = {:.3f} mm".format(result["Width"]))
    print("Coil constructive Q-factor Q ≈ {:.3f}".format(result["Q"]))
    return result


def calc_L_PCB_coil_Spiral(N=5, f=1.0, d=10.0, s=1.0, W=0.5, t=0.1):
    """calc_L_PCB_coil
    https://coil32.net/pcb-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil9.png

    Args:
        N (int): Number of turns. Defaults to 5
        f (float): Frequency f (MHz). Defaults to 1.0.
        d (float): Inside diameter (mm). Defaults to 10.0.
        s (float): Winding pitch (mm). Defaults to 1.0.
        W (float): Width of a PCB trace (mm). Defaults to 0.5.
        t (float, optional): PCB trace thickness (mm). Defaults to 0.1.

    Returns:
        _type_: result
    """

    # 0: Square; 1: Spiral; 2: Rectangular
    layoutPCB = layoutType.Spiral

    result = coil64_lib.getPCB_I(N, d, s, layoutPCB.value)
    A = d + 2 * s * N
    B = d
    result["Q"] = coil64_lib.solve_Qpcb(
        N, result["L"], A, B, W, t, s, f, layoutPCB.value)

    print("Inductance L (uH) = {:.3f} uH".format(result["L"]))
    print("Outside diameter D = {:.3f} mm".format(result["Do"]))
    print("Coil constructive Q-factor Q ≈ {:.3f}".format(result["Q"]))
    return result


# TODO: dont work!
def calc_PCB_coil_Rectangular(I=10.0, f=1.0, A=20.0, B=15.0, a=5.0, t=0.1, ratio=0.6):
    """calc_PCB_coil_Rectangular
     # TODO: dont work!
    https://coil32.net/pcb-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil8r.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        f (float): Frequency f (MHz). Defaults to 1.0.
        A (float): Outside dimension (mm). Defaults to 20.0.
        B (float): Outside dimension (mm). Defaults to 15.0.
        a (float): Inside dimension (mm). Defaults to 5.0.
        t (float, optional): PCB trace thickness (mm). Defaults to 0.1.
        ratio (float, optional): = wire width / wire pitch (mm). Defaults to 0.6.

    Returns:
        _type_: result
    """

    layoutPCB = layoutType.Rectangular

    result = coil64_lib.getPCB_RectN(I, A, B, a, t, ratio)  # TODO: dont work!
    D, d = A, a

    if ((result["Winding pitch"] != 0) and (result["Width"] != 0)):
        result["Q"] = coil64_lib.solve_Qpcb(int(result["Number turns"]), I, D, d,
                                            result["Width"], t, result["Winding pitch"],
                                            f, layoutPCB.value)
    else:
        result["Q"] = 0

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Winding pitch s = {:.3f} mm".format(result["Winding pitch"]))
    print("Width of a PCB trace W = {:.3f} mm".format(result["Width"]))
    print("Coil constructive Q-factor Q ≈ {:.3f}".format(result["Q"]))
    return result


def calc_L_PCB_coil_Rectangular(N=5, f=1.0, A=20.0, B=20.0, s=1.0, W=0.5, t=0.1):
    """calc_L_PCB_coil_Rectangular
    https://coil32.net/pcb-coil.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil8r.png

    Args:
        N (int): Number of turns. Defaults to 5.
        f (float): Frequency f (MHz). Defaults to 1.0.
        A (float): Outside dimension (mm). Defaults to 20.0.
        B (float): Outside dimension (mm). Defaults to 20.0.
        s (float): Winding pitch (mm). Defaults to 1.0.
        W (float): Width of a PCB trace (mm). Defaults to 0.5.
        t (float, optional): PCB trace thickness (mm). Defaults to 0.1.

    Returns:
        _type_: result
    """

    # 0: Square; 1: Square; 2: Rectangular
    layoutPCB = layoutType.Rectangular

    result = coil64_lib.getPCB_RectI(N, A, B, s, W, t)
    result["Q"] = coil64_lib.solve_Qpcb(
        N, result["L"], A, B, W, t, s, f, layoutPCB.value)

    print("Inductance L (uH) = {:.3f} uH".format(result["L"]))
    print("Outside diameter D = {:.3f} mm".format(result["Do"]))
    print("Coil constructive Q-factor Q ≈ {:.3f}".format(result["Q"]))
    return result


def calc_Flat_Spiral(I=10.0, Di=5.0, d=0.1, s=0.75):
    """
    https://coil32.net/foil-wound-coil-calculation.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil10.png

    Args:
        I (float): Inductance L (uH). Defaults to 10.0.
        Di (float): Inside diameter (mm). Defaults to 5.0.
        d (float): Wire diameter (mm). Defaults to 0.1.
        s (float): Gap between turns (mm). Defaults to 0.75.

    Returns:
        _type_: result
    """

    result = coil64_lib.getSpiralN(I, Di, d, s)

    print("Number of turns of the coil N = {:.3f}".format(
        result["Number turns"]))
    print("Outside diameter Do = {:.3f} mm".format(result["Do"]))
    print("Length of wire without leads lw = {:.3f} mm".format(1000 *
                                                               result["Length spiral"]))
    return result


def calc_L_Flat_Spiral(N=5, OD=25.0, ID=10.0, d=0.1):
    """calc_L_Flat_Spiral 
    https://coil32.net/foil-wound-coil-calculation.html
    https://github.com/radioacoustick/Coil64/tree/master/res/Coil10.png

    Args:
        N (int): Number of turns. Defaults to 5.
        OD (float): Outside diameter (mm). Defaults to 25.0.
        ID (float): Inside diameter (mm). Defaults to 10.0.
        d (float, optional): Wire diameter (mm). Defaults to 0.1.

    Returns:
        _type_: result
    """

    result = coil64_lib.getSpiralI(OD, ID, d, N)

    print("Inductance L = {:.3f} uH".format(result["Number turns"]))
    print("Length of wire without leads lw = {:.3f} mm".format(1000 *
                                                               result["Length spiral"]))
    return result


##########################################################################
# only plotting
##########################################################################


class data_linewidth_plot():
    def __init__(self, x, y, **kwargs):
        self.ax = kwargs.pop("ax", plt.gca())
        self.fig = self.ax.get_figure()
        self.lw_data = kwargs.pop("linewidth", 1)
        self.lw = 1
        self.fig.canvas.draw()

        self.ppd = 72./self.fig.dpi
        self.trans = self.ax.transData.transform
        self.linehandle, = self.ax.plot([], [], **kwargs)
        if "label" in kwargs:
            kwargs.pop("label")
        self.line, = self.ax.plot(x, y, **kwargs)
        self.line.set_color(self.linehandle.get_color())
        self._resize()
        self.cid = self.fig.canvas.mpl_connect('draw_event', self._resize)

    def _resize(self, event=None):
        lw = ((self.trans((1, self.lw_data))-self.trans((0, 0)))*self.ppd)[1]
        if lw != self.lw:
            self.line.set_linewidth(lw)
            self.lw = lw
            self._redraw_later()

    def _redraw_later(self):
        self.timer = self.fig.canvas.new_timer(interval=10)
        self.timer.single_shot = True
        self.timer.add_callback(lambda: self.fig.canvas.draw_idle())
        self.timer.start()


def Gen_Spule(xLength, yLength, WieVielEck, NWindungen, pitch):
    """Gen_Spule

    Args:
        xLength (_type_): xLength
        yLength (_type_): yLength
        WieVielEck (_type_): WieVielEck
        NWindungen (_type_): NWindungen
        pitch (_type_): pitch

    Returns:
        _type_: (Coilx, Coily)
    """

    if NWindungen >= (min(xLength, yLength) / pitch/2+1/2):
        NWindungen = min(xLength, yLength) / pitch/2+1/2

    NWindungen = int(NWindungen)

    if(WieVielEck == 4):
        Coilx = np.array([-0.5, 0.5, 0.5, -0.5, -0.5]) * xLength
        Coily = np.array([-0.5, -0.5, 0.5, 0.5]) * yLength
        Coily = np.append(Coily, [-0.5*yLength+pitch])

        for n in range(1, NWindungen):
            xLength = xLength - 2*pitch
            yLength = yLength - 2*pitch

            Coilx = np.append(Coilx, np.array([0.5, 0.5, -0.5, - 0.5])*xLength)
            Coily = np.append(Coily, np.array([-0.5, 0.5, 0.5])*yLength)
            Coily = np.append(Coily, [-0.5*yLength+pitch])

    elif (WieVielEck > 4):
        Coilx = np.arange(0, WieVielEck*NWindungen)
        Coily = np.arange(0, WieVielEck*NWindungen)

        if (WieVielEck > 8):
            reduction = Coilx/WieVielEck
        else:
            reduction = np.round(Coilx/WieVielEck-0.4)

        Coilx = np.sin(Coilx/WieVielEck*2 * np.pi) * \
            (xLength/2-pitch*reduction)
        Coily = np.cos(Coily/WieVielEck*2 * np.pi) * \
            (yLength/2-pitch*reduction)

    return (Coilx, Coily)


# draw_coil(num_windings=3, pitch=2, wire_width=1, outer_length=10)
def draw_coil(num_windings, pitch, wire_width, outer_length=0, inner_length=0):
    """draw_coil

    Args:
        num_windings (int): num_windings
        pitch (float): pitch
        wire_width (float): wire_width
        outer_length (int, optional): outer_length (mm). Defaults to 0.
        inner_length (int, optional): inner_length (mm). Defaults to 0.
    """

    import matplotlib.pyplot as plt

    if (outer_length == 0 and inner_length > 0):
        outer_length = inner_length + pitch * 2 * num_windings

    # Erzeugung der x- und y-Koordinatenpunkte der Spule
    x_points, y_points = Gen_Spule(
        xLength=outer_length, yLength=outer_length,
        WieVielEck=6, NWindungen=num_windings, pitch=pitch)

    # Zeichnen der Spule
    l = data_linewidth_plot(x_points, y_points, linewidth=wire_width)

    plt.axis('equal')
    plt.xlabel('Länge (mm)')
    plt.ylabel('Breite (mm)')
    plt.show()

# # TODO below
# def calc_L_Onelayer_p(N=5, f=1, D=10.0, w=1.0, t=0.0, i=0.01, p=1.5, mt=Material.Cu):
#     """One layer coil with rect wire
#     /images/res/Coil2_square.png

#     Args:
#         N (int, optional): Number of turns. Defaults to 5.
#         f (int, optional):  Frequency (MHz). Defaults to 1.
#         D (float, optional):  Former diameter (mm). Defaults to 10.0.
#         w (float, optional):  Wire width (mm). Defaults to 1.0.
#         t (float, optional):  Wire thickness (mm). Defaults to 0.0.
#         i (float, optional):  Insulation thickness (mm). Defaults to 0.01.
#         p (float, optional):  Winding pitch (mm). Defaults to 1.5.
#         mt (_type_, optional):  Material. Defaults to Material.Cu.
#     """


# # TODO below
# def calc_L_Onelayer_q(N=5, f=1, D=15.0, d=0.1, k=0.11, p=1.5, n=4, mt=Material.Cu):
#     """One layer coil on a polygonal former
#    https://coil32.net/one-layer-air-core-coil.html
#    /images/res/Coil3.png

#     Args:
#         N (int, optional): Number of turns. Defaults to 5.
#         f (int, optional): Frequency (MHz). Defaults to 1.
#         D (float, optional): Former diameter (mm). Defaults to 15.0.
#         d (float, optional): Wire diameter (mm). Defaults to 0.1.
#         k (float, optional): Wire diameter with insulation (mm). Defaults to 0.11.
#         p (float, optional): Winding pitch (mm). Defaults to 1.5.
#         n (int, optional): Number of sides of the former. Defaults to 4.
#         mt (_type_, optional): Material. Defaults to Material.Cu.
#     """


# # TODO below
# def calc_L_Multilayer(N=5, D=15.0, l=8.0, d=0.1, k=0.11):
#     """Multilayer coil
#    https://coil32.net/multi-layer-coil.html
#    /images/res/Coil4.png

#     Args:
#         N (int, optional): Number of turns. Defaults to 5.
#         D (float, optional): Former diameter (mm). Defaults to 15.0.
#         l (float, optional): Winding length (mm). Defaults to 8.0.
#         d (float, optional): Wire diameter (mm). Defaults to 0.1.
#         k (float, optional): Wire diameter with insulation (mm). Defaults to 0.11.
#     """
#     # c = 0.0  # Thickness of the coil (mm)


# # TODO below
# def calc_L_Multilayer_p(N=5, D=15.0, l=8.0, c=1.00, d=0.1, k=0.11, g=0.1, Ng=5):
#     """Multilayer coil with insulated pads
#     # https://coil32.net/multi-layer-coil.html
#     # /images/res/Coil4-0.png
#     Args:
#         N (int, optional): Number of turns. Defaults to 5.
#         D (float, optional): Former diameter (mm). Defaults to 15.0.
#         l (float, optional): Winding length (mm). Defaults to 8.0.
#         c (float, optional): Thickness of the coil (mm). Defaults to 1.00.
#         d (float, optional): Wire diameter (mm). Defaults to 0.1.
#         k (float, optional): Wire diameter with insulation (mm). Defaults to 0.11.
#         g (float, optional): Insulation thickness (mm). Defaults to 0.1.
#         Ng (int, optional): Layers number beetween insulating pads. Defaults to 5.
#     """


# # TODO below
# def calc_L_Multilayer_r(N=5, D=15.0, b=2.0, l=8.0, c=1.00, d=0.1, k=0.11):
#     """Multilayer coil on a rectangular former
#     # https://coil32.net/multilayer-rectangular.html
#     # /images/res/Coil4_square.png

#     Args:
#         N (int, optional): Number of turns. Defaults to 5.
#         D (float, optional): Former diameter (mm). Defaults to 15.0.
#         b (float, optional): Former height (mm). Defaults to 2.0.
#         l (float, optional): Winding length (mm). Defaults to 8.0.
#         c (float, optional): Thickness of the coil (mm). Defaults to 1.00.
#         d (float, optional): Wire diameter (mm). Defaults to 0.1.
#         k (float, optional): Wire diameter with insulation (mm). Defaults to 0.11.
#     """
