

from matplotlib import patches
import matplotlib.pyplot as plt

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


def calc_Multilayer():  # Multilayer coil
    # https://coil32.net/multi-layer-coil.html
    # /images/res/Coil4.png

    I = 20  # Inductance L
    D = 8  # Former diameter
    lk = 5  # Winding length
    d = 0.1  # Wire diameter
    k = 0.109  # Wire diameter with insulation

    result = coil64_lib.getMultiLayerN(I, D, d, k, lk, 0, -1)

    print("Number of turns of the coil N = {}".format(
        result['Number turns']))  # TODO: can be wrong!
    print("Thickness of the coil c = {} mm".format(result['Thickness']))
    print("Length of wire without leads lw = {} m".format(result['Length']))
    print("DC resistance of the coil Rdc = {} Ohm".format(result['R_DC']))
    print("Number of layers Nl = {}".format(result['Number layers']))


def calc_Multilayer_p():  # Multilayer coil with insulated pads
    # https://coil32.net/multi-layer-coil.html
    # /images/res/Coil4-0.png

    I = 10  # Inductance L
    D = 15  # Former diameter
    lk = 8  # Winding length
    d = 0.1  # Wire diameter
    k = 0.11  # Wire diameter with insulation
    g = 0.1  # Insulation thickness
    Ng = 5  # Layers number beetween insulating pads

    result = coil64_lib.getMultiLayerN(I, D, d, k, lk, g, Ng)

    print("Number of turns of the coil N = {}".format(result["Number turns"]))
    print("Thickness of the coil c = {} mm".format(result["Thickness"]))
    print("Length of wire without leads lw = {} m".format(result["Length"]))
    print("DC resistance of the coil Rdc = {} Ohm".format(result["R_DC"]))
    print("Number of layers Nl = {}".format(result["Number layers"]))
    print("Number of interlayers Ng = {}".format(result["Ng"]))


def calc_Multilayer_r():  # Multilayer coil on a rectangular former
    # https://coil32.net/multilayer-rectangular.html
    # /images/res/Coil4_square.png

    I = 10  # Inductance L
    a = 5  # Former width
    b = 2  # Former height
    l = 8  # Winding length
    d = 0.1  # Wire diameter
    k = 0.11  # Wire diameter with insulation

    result = coil64_lib.getMultiLayerN_rectFormer(I, a, b, l, d, k)

    print("Number of turns of the coil N = {}".format(result["Number turns"]))
    print("Number of layers Nl = {}".format(result["Number layers"]))
    print("Thickness of the coil c = {} mm".format(result["thickness"]))
    print("Length of wire without leads lw = {} m".format(
        result["Length wire"]))
    print("DC resistance of the coil Rdc = {} Ohm".format(result["Rdc"]))


def calc_Multilayer_f():  # Multilayer foil-wound coil
    # https://coil32.net/foil-wound-coil-calculation.html
    # /images/res/Coil11.png

    I = 10  # Inductance L
    D = 15  # Former diameter
    w = 12  # Foil width
    t = 0.1  # Foil thickness
    g = 0.01  # Insulation thickness
    ins = g

    result = coil64_lib.getMultilayerN_Foil(D, w, t, ins, I)

    print("Number of turns of the coil N = {}".format(result["Number turns"]))
    print("Outside diameter Do = {} mm".format(result["Do"]))
    print("Length of the foil lf = {} m".format(result["Length spiral"]))
    print("DC resistance of the coil Rdc = {} Ohm (Copper)".format(
        result["Rdcc"]))
    print("DC resistance of the coil Rdc = {} Ohm (Aluminum)".format(
        result["Rdca"]))


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

    result = coil64_lib.getFerriteN(I, OD, ID, h, d, mu, C)

    print("Number of turns of the coil N = {}".format(result["Number turns"]))
    print("Length of wire without leads lw = {} m".format(
        result["Length wire"]))
    print("AL = {} nH/N2".format(result["Al"]))


def calc_PCB_coil():  # only Square and Spiral are working!
    # https://coil32.net/pcb-coil.html

    I = 10.0  # Inductance L
    f = 1.0  # Frequency f

    # 0: square; 1: spiral; 2: rectangular
    layoutPCB = layoutType.Square

    if layoutPCB == layoutType.Square:
        # /images/res/Coil8.png
        D = 20.0  # Outside diameter
        d = 5.0  # Inside diameter
        t = 0.1  # PCB trace thickness
    elif layoutPCB == layoutType.Spiral:
        # /images/res/Coil9.png
        D = 20.0  # Outside diameter
        d = 5.0  # Inside diameter
        t = 0.1  # PCB trace thickness
    elif layoutPCB == layoutType.Rectangular:
        # /images/res/Coil8r.png
        A = 20.0  # Outside dimension
        B = 15.0  # Outside dimension
        a = 5.0  # Inside dimension
        t = 0.1  # PCB trace thickness

    ratio = 0.6  # = wire width / wire pitch

    if (layoutPCB != layoutType.Rectangular):
        result = coil64_lib.calc_getPCB_N(I, D, d, ratio, layoutPCB.value)
    else:
        result = coil64_lib.getPCB_RectN(
            I, A, B, a, t, ratio)  # TODO: dont work!
        D, d = A, a

    if ((result["Winding pitch"] != 0) and (result["Width"] != 0)):
        result["Q"] = coil64_lib.solve_Qpcb(int(result["Number turns"]), I, D, d,
                                            result["Width"], t, result["Winding pitch"],
                                            f, layoutPCB.value)
    else:
        result["Q"] = 0

    print("Number of turns of the coil N = {}".format(result["Number turns"]))
    print("Winding pitch s = {}mm".format(result["Winding pitch"]))
    print("Width of a PCB trace W = {} mm".format(result["Width"]))
    print("Coil constructive Q-factor Q ≈ {}".format(result["Q"]))


def calc_Flat_Spiral():
    # https://coil32.net/foil-wound-coil-calculation.html
    # /images/res/Coil10.png

    I = 10  # Inductance L
    Di = 5  # Inside diameter
    d = 0.1  # Wire diameter
    s = 0.75  # Gap between turns

    result = coil64_lib.getSpiralN(I, Di, d, s)

    print("Number of turns of the coil N = {}".format(result["Number turns"]))
    print("Outside diameter Do = {} mm".format(result["Do"]))
    print("Length of wire without leads lw = {} m".format(
        result["Length spiral"]))


# calc_Multilayer()
# calc_Multilayer_p()
# calc_Multilayer_r()
# calc_Multilayer_f()
# calc_FerrToroid()
# calc_PCB_coil()
# calc_Flat_Spiral()

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


def draw_coil(num_windings, pitch, wire_width, outer_length=0, inner_length=0):

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


# draw_coil(num_windings=3, pitch=2, wire_width=1, outer_length=10)
