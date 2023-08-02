import lib_python.resolves as resolves
import lib_python.resolve_q as resolve_q
import lib_python.resolve_srf_cs as resolve_srf_cs

result = resolves._CoilResult()


# is working # _Onelayer_cw  # One layer close-winding coil
# is working # _Onelayer     # One layer coil with round wire
# not working # _Onelayer_p   # One layer coil with rect wire
# not working # _Onelayer_q   # One layer coil on a polygonal former
# not working # _Multilayer   # Multilayer coil
# not working # _Multilayer_p # Multilayer coil with insulated pads
# not working # _Multilayer_r # Multilayer coil on a rectangular former
# is working # _Multilayer_f # Multilayer foil-wound coil
# is  working # _FerrToroid   # Ferrite toroid coil
# ok  working # _PCB_coil     # PCB flat coil
# is  working # _Flat_Spiral  # Tesla flat spiral coil

def calc_L_Onelayer_cw():  # One layer close-winding coil
    # /images/res/Coil1.png

    N = 0  # Number of turns
    f = 0  # Frequency
    D = 0  # Former diameterd
    d = 0  # Wire diameterd
    k = 0  # Wire diameter with insulationd

    # same like calc_L_Onelayer()


def calc_L_Onelayer():  # One layer coil with round wire
    # /images/res/Coil2.png

    N = 15  # Number of turns
    f = 1  # Frequency
    D = 15  # Former diameterd
    d = 0.1  # Wire diameterd
    k = 0.11  # Wire diameter with insulationd # TODO to check why unused
    p = 1.5  # Winding pitch

    mt = resolve_q.Material.Cu
    accuracy = 1

    L_uH, lw = resolves.getOneLayerI_withRoundWire(D, d, p, N, accuracy)
    result.thd = resolve_srf_cs.find_Cs(p, D, p * N)  # self-capacitance
    result.fourth = resolve_srf_cs.findSRF(
        p * N, D, lw)  # self-resonance frequency
    QualitFactor = resolve_q.solve_Qr(
        L_uH, D, p, d, f, N, result.thd, mt, result)  # Q-factor

    print("Inductance L = {} uH".format(L_uH))
    print("Length of wire without leads lw = {} m".format(lw))
    print("Self capacitance Cs = {} pF".format(result.thd))
    print("Coil self-resonance frequency Fsr = {} MHz".format(result.fourth))
    print("Coil constructive Q-factor Q = {}".format(QualitFactor))
    print("Loss resistance ESR = {} Ohm".format(result.seven))


def calc_L_Onelayer_p():  # One layer coil with rect wire
    # /images/res/Coil2_square.png

    N = 0  # Number of turns
    f = 0  # Frequency
    D = 0  # Former diameterd
    w = 0  # Wire widthd
    t = 0  # Wire thicknessd
    i = 0  # Insulation thickness
    p = 0  # Winding pitch


def calc_L_Onelayer_q():  # One layer coil on a polygonal former
    # https://coil32.net/one-layer-air-core-coil.html
    # /images/res/Coil3.png

    N = 0  # Number of turns
    f = 0  # Frequency
    D = 0  # Former diameterd
    d = 0  # Wire diameterd
    k = 0  # Wire diameter with insulationd
    p = 0  # Winding pitch
    n = 0  # Number of sides of the former


def calc_L_Multilayer():  # Multilayer coil
    # https://coil32.net/multi-layer-coil.html
    # /images/res/Coil4.png

    N = 0  # Number of turns
    D = 0  # Former diameterd
    l = 0  # Winding lengthd
    c = 0  # Thickness of the coild
    d = 0  # Wire diameterd
    k = 0  # Wire diameter with insulationd


def calc_L_Multilayer_p():  # Multilayer coil with insulated pads
    # https://coil32.net/multi-layer-coil.html
    # /images/res/Coil4-0.png

    N = 0  # Number of turns
    D = 0  # Former diameterd
    l = 0  # Winding lengthd
    c = 0  # Thickness of the coild
    d = 0  # Wire diameterd
    k = 0  # Wire diameter with insulationd
    g = 0  # Insulation thickness
    Ng = 0  # Layers number beetween insulating padsd


def calc_L_Multilayer_r():  # Multilayer coil on a rectangular former
    # https://coil32.net/multilayer-rectangular.html
    # /images/res/Coil4_square.png

    N = 0  # Number of turns
    a = 0  # Former widthd
    b = 0  # Former heightd
    l = 0  # Winding lengthd
    c = 0  # Thickness of the coild
    d = 0  # Wire diameterd
    k = 0  # Wire diameter with insulationd


def calc_L_Multilayer_f():  # Multilayer foil-wound coil
    # https://coil32.net/foil-wound-coil-calculation.html
    # /images/res/Coil11.png

    N = 10  # Number of turns
    D = 15  # Former diameterd
    w = 10  # Foil widthd
    t = 0.1  # Foil thicknessd
    g = 0.01  # Insulation thickness

    ins = g

    resolves.getMultilayerI_Foil(D, w, t, ins, N, result)

    print(result)

    print("Inductance L = {} uH".format(result.N))
    print("Length of the foil lf = {} m".format(result.sec))
    print("Outside diameter Do = {} mm".format(result.thd))
    print("DC resistance of the coil Rdc = {} Ohm (Copper)".format(result.fourth))
    print("DC resistance of the coil Rdc = {} Ohm (Aluminum)".format(result.five))


def calc_L_FerrToroid():
    # https://coil32.net/ferrite-toroid-core.html
    # /images/res/T-core.png

    N = 5  # Number of turns
    OD = 20  # Outside diameter
    ID = 10  # Inside diameter
    h = 5  # Core heightd
    mu = 10  # Init magnetic permeability
    C = 0  # Chamferd

    L_uH = resolves.getFerriteI(N, OD, ID, h, mu, C, result)

    print("Inductance L = {} uH".format(L_uH))
    print("A_L = {} nH/N^2".format(result.thd))


def calc_L_PCB_coil():
    # https://coil32.net/pcb-coil.html

    N = 5  # Number of turns
    f = 1  # Frequency

    # 0: square; 1: spiral; 2: rectangular
    layoutPCB = resolve_q.layoutType.Square

    if layoutPCB == resolve_q.layoutType.Square:
        # /images/res/Coil8.png
        d = 10  # Inside diameter
        s = 1  # Winding pitch
        W = 0.5  # Width of a PCB trace
        t = 1  # PCB trace thickness
    elif layoutPCB == resolve_q.layoutType.Spiral:
        # /images/res/Coil8.png
        d = 10  # Inside diameter
        s = 1  # Winding pitch
        W = 0.5  # Width of a PCB trace
        t = 1  # PCB trace thickness
    elif layoutPCB == resolve_q.layoutType.Rectangular:
        # /images/res/Coil8r.png
        A = 20  # Outside dimension
        B = 20  # Outside dimension
        s = 1  # Winding pitch
        W = 0.5  # Width of a PCB trace
        t = 0.1  # PCB trace thickness

    if (layoutPCB != resolve_q.layoutType.Rectangular):
        L_uH = resolves.getPCB_I(N, d, s, layoutPCB.value, result)
        A = d + 2 * s * N  # TODO: incorrect
        B = d
    else:
        L_uH = resolves.getPCB_RectI(N, A, B, s, W, t, result)

    QualitFactor = resolve_q.solve_Qpcb(N, L_uH, A, B, W, t, s, f, layoutPCB)

    print("Inductance L = {} uH".format(L_uH))
    print("Outside diameter D = {} mm".format(result.five))
    print("Coil constructive Q-factor Q ≈ {}".format(QualitFactor))


def calc_L_Flat_Spiral():
    # https://coil32.net/foil-wound-coil-calculation.html
    # /images/res/Coil10.png

    N = 0  # Number of turns
    OD = 25  # Outside diameter
    ID = 10  # Inside diameter
    d = 0.1  # Wire diameter

    resolves.getSpiralI(OD, ID, d, N, result)

    print("Inductance L = {} uH".format(result.N))
    print("Length of wire without leads lw = {} m".format(result.sec))
