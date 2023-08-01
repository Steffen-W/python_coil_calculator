import lib_python.resolves as resolves
import lib_python.resolve_q as resolve_q
import lib_python.resolve_srf_cs as resolve_srf_cs

result = resolves._CoilResult()


# not working # _Onelayer_cw  # One layer close-winding coil
# not working # _Onelayer     # One layer coil with round wire
# is  working # _Onelayer_p   # One layer coil with rect wire
# is  working # _Onelayer_q   # One layer coil on a polygonal former
# is  working # _Multilayer   # Multilayer coil
# is  working # _Multilayer_p # Multilayer coil with insulated pads
# is  working # _Multilayer_r # Multilayer coil on a rectangular former
# not working # _Multilayer_f # Multilayer foil-wound coil
# is  working # _FerrToroid   # Ferrite toroid coil
# ok  working # _PCB_coil     # PCB flat coil
# is  working # _Flat_Spiral  # Tesla flat spiral coil


# ERROR: currently unusable
def calc_Onelayer_cw():  # One layer close-winding coil
    # /images/res/Coil1.png
    I = 5  # Inductance L
    f = 1  # Frequency f

    D = 100  # Former diameter
    d = 0.1  # Wire diameter
    k = 0.11  # Wire diameter with insulation

    RoundWire = False
    l = 5  # Winding length

    mt = resolve_q.Material.Cu
    accuracy = 1

    arg1, arg2, arg3, arg4, arg5 = D, d, k, I, f
    Dk = 0
    dw = 0
    if RoundWire:  # when winding length option is activated, k=0
        # result.N, result.sec = resolves.getOneLayerN_withRoundWire(
        #     arg1, arg2, arg3, arg4, accuracy)  # number of turns
        Dk = arg1
        dw = arg2
    else:
        result.N = resolves.getOneLayerN_byWindingLength(
            arg1, arg2, arg4, result, accuracy)
        dw = result.five
        arg3 = resolves.odCalc(dw)
        Dk = arg1 + arg3

    result.thd = resolve_srf_cs.find_Cs(
        arg3, Dk, arg3 * result.N)  # self-capacitance
    result.six = resolve_q.solve_Qr(
        arg4, Dk, arg3, dw, arg5, result.N, result.thd, mt, result)  # Q-factor
    # result.fourth = resolve_srf_cs.findSRF(
    #     arg3 * result.N, Dk, result.sec)  # self-resonance frequency

    print(result)
    print(Dk, dw)

    print("Number of turns of the coil {}".format(result.N))
    # print("Length of wire without leads lw = {} m".format(result.six))
    # print("Length of winding l = {} mm".format(result.))
    # print("DC resistance of the coil Rdc = {} Ohm".format(result.))
    print("Reactance of the coil X = {} Ohm".format(result.thd))

    print("Self capacitance Cs = {} pF".format(result.fourth))
    print("Coil self-resonance frequency Fsr = {} MHz".format(result.fourth))
    print("Coil constructive Q-factor = {}".format(result.six))
    print("Loss resistance ESR = {} Ohm".format(result.seven))


# ERROR: currently unusable
def calc_Onelayer():  # One layer coil with round wire
    # /images/res/Coil2.png
    I = 5  # Inductance L
    f = 1  # Frequency f
    D = 0  # Former diameter
    d = 0  # Wire diameter
    k = 0  # Wire diameter with insulation
    p = 0  # Winding pitch

    mt = resolve_q.Material.Cu
    # calc_Onelayer_cw()


def calc_Onelayer_p():  # One layer coil with rect wire
    # /images/res/Coil2_square.png
    I = 10  # Inductance L
    f = 1  # Frequency f

    D = 10  # Former diameter
    w = 1  # Wire width
    t = 0.1  # Wire thickness
    i = 0.1  # Insulation thickness
    p = 1.5  # Winding pitch

    mt = resolve_q.Material.Cu
    accuracy = 1

    Dk = D
    t = i

    result.N, result.sec = resolves.getOneLayerN_withRectWire(
        D, w, t, p, I, accuracy)  # number of turns
    result.thd = resolve_srf_cs.find_Cs(
        p, D, p * result.N)  # self-capacitance
    result.six = resolve_q.solve_Qc(
        I, D, p, w, t, f, result.N, result.thd, mt, result)  # Q-factor
    result.fourth = resolve_srf_cs.findSRF(
        p * result.N, D, result.sec)  # self-resonance frequency

    print("Number of turns of the coil N = {}".format(result.N))
    print("Length of wire without leads lw = {} m".format(result.sec))
    # print("Length of winding l = 115.43 mm".format(result.))
    # print("DC resistance of the coil Rdc = 0.418 Ohm".format(result.))
    # print("Reactance of the coil X = 31.416 Ohm".format(result.))

    print("Self capacitance Cs = 1.418 pF".format(result.thd))
    print("Coil self-resonance frequency Fsr = {} MHz".format(result.fourth))
    print("Coil constructive Q-factor Q = {}".format(result.six))
    print("Loss resistance ESR = {} Ohm".format(result.seven))


def calc_Onelayer_q():  # One layer coil on a polygonal former
    # https://coil32.net/one-layer-air-core-coil.html
    # /images/res/Coil3.png

    I = 10  # Inductance L
    f = 1  # Frequency f
    D = 10  # Former diameter
    d = 0.1  # Wire diameter
    k = 0.11  # Wire diameter with insulation
    p = 1.5  # Winding pitch
    n = 6  # Number of sides of the former

    mt = resolve_q.Material.Cu
    accuracy = 1

    # number of turns
    result.N = resolves.getOneLayerN_Poligonal(I, D, d, p, n, result, accuracy)
    result.fourth = resolve_srf_cs.find_Cs(
        p, result.seven, p * result.N)    # self-capacitance
    # self-resonance frequency    Q-factor
    result.five = resolve_srf_cs.findSRF(
        p * result.N, result.seven, result.thd)
    result.six = resolve_q.solve_Qr(
        I, result.seven, p, d, f, result.N, result.fourth, mt, result)

    print("Number of turns of the coil N = {}".format(result.N))
    print("Length of wire without leads lw = {} m".format(result.thd))
    print("Length of winding l = {} mm".format(result.sec))
    # print("DC resistance of the coil Rdc = {} Ohm".format(result.))
    # print("Reactance of the coil X = {} Ohm".format(result.))
    print("Self capacitance Cs = {} pF".format(result.fourth))
    print("Coil self-resonance frequency Fsr = {} MHz".format(result.five))
    print("Coil constructive Q-factor Q = {}".format(result.six))
    print("Loss resistance ESR = {} Ohm".format(result.seven))


def calc_Multilayer():  # Multilayer coil
    # https://coil32.net/multi-layer-coil.html
    # /images/res/Coil4.png

    I = 10  # Inductance L
    D = 10  # Former diameter
    lk = 5  # Winding length
    d = 0.1  # Wire diameter
    k = 0.109  # Wire diameter with insulation

    resolves.getMultiLayerN(I, D, d, k, lk, 0, -1, result)

    print("Number of turns of the coil N = {}".format(result.six))
    print("Thickness of the coil c = {} mm".format(result.fourth))
    print("Length of wire without leads lw = {} m".format(result.sec))
    print("DC resistance of the coil Rdc = {} Ohm".format(result.N))
    print("Number of layers Nl = {}".format(result.thd))


def calc_Multilayer_p():  # Multilayer coil with insulated pads
    # https://coil32.net/multi-layer-coil.html
    # /images/res/Coil4-0.png

    I = 10  # Inductance L
    D = 10  # Former diameter
    lk = 5  # Winding length
    d = 0.1  # Wire diameter
    k = 0.11  # Wire diameter with insulation
    g = 0.1  # Insulation thickness
    Ng = 5  # Layers number beetween insulating pads

    resolves.getMultiLayerN(I, D, d, k, lk, g, Ng, result)

    print("Number of turns of the coil N = {}".format(result.six))
    print("Thickness of the coil c = {} mm".format(result.fourth))
    print("Length of wire without leads lw = {} m".format(result.sec))
    print("DC resistance of the coil Rdc = {} Ohm".format(result.N))
    print("Number of layers Nl = {}".format(result.thd))
    print("Number of interlayers Ng = {}".format(result.five))


def calc_Multilayer_r():  # Multilayer coil on a rectangular former
    # https://coil32.net/multilayer-rectangular.html
    # /images/res/Coil4_square.png

    I = 10  # Inductance L
    a = 1  # Former width
    b = 2  # Former height
    l = 5  # Winding length
    d = 0.1  # Wire diameter
    k = 0.11  # Wire diameter with insulation

    resolves.getMultiLayerN_rectFormer(I, a, b, l, d, k, result)

    print("Number of turns of the coil N = {}".format(result.N))
    print("Number of layers Nl = {}".format(result.sec))
    print("Thickness of the coil c = {} mm".format(result.five))
    print("Length of wire without leads lw = {} m".format(result.thd))
    print("DC resistance of the coil Rdc = 1.332 Ohm".format(result.fourth))


# dont work correct !
def calc_Multilayer_f():  # Multilayer foil-wound coil
    # https://coil32.net/foil-wound-coil-calculation.html
    # /images/res/Coil11.png

    I = 10  # Inductance L
    D = 100  # Former diameter
    w = 10  # Foil width
    t = 0.01  # Foil thickness
    g = 0.1  # Insulation thickness

    # ins = g

    # resolves.getMultilayerN_Foil(D, w, t, ins, I, result)

    # Number of turns of the coil N = {}
    # Outside diameter Do = {} mm
    # Length of the foil lf = {} m
    # DC resistance of the coil Rdc = {} Ohm (Copper)
    # DC resistance of the coil Rdc = {} Ohm (Aluminum)

    print(result)  # dont work correct !


def calc_FerrToroid():
    # https://coil32.net/ferrite-toroid-core.html
    # /images/res/T-core.png

    I = 10  # Inductance L
    OD = 30  # Outside diameter
    ID = 10  # Inside diameter
    h = 5  # Core height
    d = 0.1  # Wire diameter
    mu = 10  # Init magnetic permeability
    C = 0  # Chamfer

    resolves.getFerriteN(I, OD, ID, h, d, mu, C, result)

    print("Number of turns of the coil N = {}".format(result.N))
    print("Length of wire without leads lw = {} m".format(result.sec))
    print("AL = {} nH/N2".format(result.thd))


# only Square and Spiral are working!
def calc_PCB_coil():
    # https://coil32.net/pcb-coil.html

    I = 6  # Inductance L
    f = 1  # Frequency f

    # 0: square; 1: spiral; 2: rectangular
    layoutPCB = resolve_q.layoutType.Square

    if layoutPCB == resolve_q.layoutType.Square:
        # /images/res/Coil8.png
        D = 20  # Outside diameter
        d = 5  # Inside diameter
        t = 0.1  # PCB trace thickness
    elif layoutPCB == resolve_q.layoutType.Spiral:
        # /images/res/Coil9.png
        D = 20  # Outside diameter
        d = 5  # Inside diameter
        t = 0.1  # PCB trace thickness
    elif layoutPCB == resolve_q.layoutType.Rectangular:
        # /images/res/Coil8r.png
        A = 20  # Outside dimension
        B = 15  # Outside dimension
        a = 5  # Inside dimension
        t = 0.1  # PCB trace thickness

    ratio = 0.6  # = wire width / wire pitch

    if (layoutPCB != resolve_q.layoutType.Rectangular):
        resolves.getPCB_N(I, D, d, ratio, layoutPCB.value, result)
    else:
        # resolves.getPCB_RectN(I, A, B, a, t, ratio, result)
        D, d = A, a

    if ((result.sec != 0) and (result.thd != 0)):
        result.fourth = resolve_q.solve_Qpcb(
            result.N, I, D, d, result.thd, t, result.sec, f, layoutPCB)

    print("Number of turns of the coil N = {}".format(result.N))
    print("Winding pitch s = {}mm".format(result.sec))
    print("Width of a PCB trace W = {} mm".format(result.thd))
    print("Coil constructive Q-factor Q ≈ {}".format(result.five))


def calc_Flat_Spiral():
    # https://coil32.net/foil-wound-coil-calculation.html
    # /images/res/Coil10.png

    I = 0.311  # Inductance L
    Di = 22.5  # Inside diameter
    d = 0.1  # Wire diameter
    s = 0.75  # Gap between turns

    resolves.getSpiralN(I, Di, d, s, result)

    print("Number of turns of the coil N = {}".format(result.N))
    print("Outside diameter Do = {} mm".format(result.thd))
    print("Length of wire without leads lw = {} m".format(result.sec))
