from enum import Enum
import math
import lib_python.resolve_q as resolve_q

# based on https://coil32.net/pcb-coil.html

µ0 = 4 * math.pi * 1e-7


class CoilType(Enum):
    Square = 0
    Hexagonal = 1
    Octagonal = 2
    Circle = 3


def calc_avg_D(D_outer: float, D_inner: float) -> float:
    D_avg = (D_outer + D_inner) / 2  # average diameter of coil
    phi = (D_outer-D_inner) / (D_outer+D_inner)  # fill factor
    return D_avg, phi


def induc_planar_spiral3(type: CoilType, N_turns: int, D_outer: float, D_inner: float) -> float:

    if type == CoilType.Square:
        c1, c2, c3, c4 = 1.27, 2.07, 0.18, 0.13
    elif type == CoilType.Hexagonal:
        c1, c2, c3, c4 = 1.09, 2.23, 0.00, 0.17
    elif type == CoilType.Octagonal:
        c1, c2, c3, c4 = 1.07, 2.29, 0.00, 0.19
    elif type == CoilType.Circle:
        c1, c2, c3, c4 = 1.00, 2.46, 0.00, 0.20

    D_avg, phi = calc_avg_D(D_outer, D_inner)

    part1 = µ0 * N_turns**2 * (D_avg*1e-3) * c1 / 2
    part2 = math.log(c2/phi)+c3*phi+c4*phi**2
    L_uH = part1 * part2

    return L_uH * 1e6


def induc_planar_spiral2(type: str, N_turns: int, D_outer: float, D_inner: float, width: float, gap: float) -> float:

    if type == CoilType.Square:
        beta, a1, a2, a3, a4, a5 = 1.62*1e-3, -1.21, -0.147, 2.40, 1.78, -0.03
    elif type == CoilType.Hexagonal:
        beta, a1, a2, a3, a4, a5 = 1.28*1e-3, -1.24, -0.174, 2.47, 1.77, -0.049
    elif type == CoilType.Octagonal:
        beta, a1, a2, a3, a4, a5 = 1.33*1e-3, -1.21, -0.163, 2.43, 1.75, -0.049

    D_avg, phi = calc_avg_D(D_outer, D_inner)

    L_uH = 1e3 * beta * (D_outer*1e-3)**a1 * (width*1e-3)**a2 * \
        (D_avg*1e-3)**a3 * N_turns**a4 * (gap*1e-3)**a5

    return L_uH


# induc_planar_spiral is not correct
def induc_planar_spiral(type: str, N_turns: int, D_outer: float, D_inner: float) -> float:

    if type == CoilType.Square:
        K1, K2 = 2.34, 2.75
    elif type == CoilType.Hexagonal:
        K1, K2 = 2.33, 3.82
    elif type == CoilType.Octagonal:
        K1, K2 = 2.25, 3.55
    elif type == CoilType.Circle:
        K1, K2 = 2.23, 3.45

    D_avg, phi = calc_avg_D(D_outer, D_inner)

    L_uH = K1 * µ0 * N_turns**2 * (D_avg*1e-3) / (1+K2*phi) * 1e6

    return L_uH


N_turns = 6
D_outer = 30
width = 0.25
gap = 0.5
pitch = gap + width
D_inner = D_outer - pitch*(N_turns-1) * 2
thickness = 0.1

freq_MHz = 1

print("N_turns  ", N_turns)
print("D_outer  ", D_outer)
print("D_inner  ", D_inner)
print("width    ", width)
print("gap      ", gap)
print("pitch    ", pitch)
print("thickness", thickness)
print("freq_MHz ", freq_MHz)

L_uH = induc_planar_spiral(CoilType.Octagonal, N_turns, D_outer, D_inner)
print("L_uH     ", L_uH)

L_uH = induc_planar_spiral2(
    CoilType.Octagonal, N_turns, D_outer, D_inner, width, gap)
print("L_uH     ", L_uH)

L_uH = induc_planar_spiral3(CoilType.Octagonal, N_turns, D_outer, D_inner)
print("L_uH     ", L_uH)

mt = resolve_q.Material.Cu

# /// Q-FACTOR OF THE PCB SPIRAL COIL
q_factor = resolve_q.solve_Qpcb(N=N_turns, I=L_uH, D=D_outer, d=D_inner, W=width,
                                t=thickness, s=pitch,  f=freq_MHz, layout=resolve_q.layoutType.Spiral)

print("q_factor ", q_factor)
