import resolves


result = resolves._CoilResult()


def calc_Multilayer():
    # /images/res/Coil4.png

    # arg: I, D, d, k, lk
    I = 0.311  # Inductance L
    D = 10  # Former diameter
    lk = 0.486  # Winding length
    d = 0.1  # Wire diameter
    k = 0.109  # Wire diameter with insulation

    resolves.getMultiLayerN(I, D, d, k, lk, 0, -1, result)

    print("Number of turns of the coil N = {}".format(result.six))
    print("Thickness of the coil c = {} mm".format(result.fourth))
    # print("Dimensions of inductor: 0.486x11x11 mm".format(result.))
    print("Length of wire without leads lw = {} m".format(result.sec))
    print("DC resistance of the coil Rdc = {} Ohm".format(result.N))
    # print("Weight of wire m = 0.016 g".format(result.))
    print("Number of layers Nl = {}".format(result.thd))


def calc_Multilayer_p():
    # /images/res/Coil4-0.png

    # arg: I, D, d, k, lk
    I = 0.311  # Inductance L
    D = 10  # Former diameter
    lk = 0.486  # Winding length
    d = 0.1  # Wire diameter
    k = 0.109  # Wire diameter with insulation

    g = 0  # Insulation thickness
    Ng = 0  # Layers number beetween insulating pads


def calc_Multilayer_r():
    # /images/res/Coil4_square.png
    I = 0.311  # Inductance L
    a = 0  # Former width
    b = 0  # Former height
    l = 0  # Winding length
    d = 0  # Wire diameter
    k = 0  # Wire diameter with insulation


def calc_Onelayer_q():
    # /images/res/Coil2.png
    # /images/res/Coil2_square.png
    # /images/res/Coil3.png
    D = 0  # Former diameter"
    w = 0  # Wire width
    t = 0  # Wire thickness
    i = 0  # Insulation thickness
    p = 0  # Winding pitch


def calc_Multilayer_f():
    # /images/res/Coil11.png
    I = 10  # Inductance L
    D = 100  # Former diameter
    w = 10  # Foil width
    t = 0.01  # Foil thickness
    g = 0.1  # Insulation thickness

    resolves.getMultilayerN_Foil(D, w, t, g, I, result)

    print(result)  # dont work correct !


def calc_FerrToroid():
    # /images/res/T-core.png
    I = 0.311  # Inductance L
    OD = 25.9  # Outside diameter
    ID = 22.5  # Inside diameter
    h = 5  # Core height
    d = 0.1  # Wire diameter
    mu = 10  # Init magnetic permeability
    C = 0  # Chamfer

    resolves.getFerriteN(I, OD, ID, h, d, mu, C, result)

    print(result)  # dont work correct !


def calc_Flat_Spiral():
    # /images/res/Coil10.png
    I = 0.311  # Inductance L
    Di = 22.5  # Inside diameter
    d = 0.1  # Wire diameter
    s = 0.75  # Gap between turns

    resolves.getSpiralN(I, Di, d, s, result)

    print("Number of turns of the coil N = {}".format(result.N))
    print("Outside diameter Do = {} mm".format(result.thd))
    print("Length of wire without leads lw = {} m".format(result.sec))
