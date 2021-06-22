import numpy as np
from scipy import constants as const

def Qf2σ(Qs, Q0, fs, f0, Vs, Vc, A = 1, B = 2, C = 0.07, delta = 1):
    p = {}
    p["ε'p"] = (1 / A) * C * (Vc/Vs) * ((f0 - fs)/fs) + 1
    p["ε'b"] = ((p["ε'p"]**(1/3) - 1)/delta + 1)**3
    p['ε"p'] = (1 / B) * C * (Vc/Vs) * (1/Qs - 1/Q0)
    p['ε"b'] = (p['ε"p']/delta) * (p["ε'b"]/p["ε'p"])**(2/3)
    p["σ"] = p['ε"b'] * const.epsilon_0 * fs * 2 * np.pi
    return p

def prettifySample(name):
    name = name.replace("a-", "α-").replace("b-", "β-")
    if "0." in name:
        name = name.replace("0.80", "$_{0.80}$").replace("0.70", "$_{0.70}$").replace("0.60", "$_{0.60}$")
        name = name.replace("0.50", "$_{0.50}$").replace("0.40", "$_{0.40}$").replace("0.60", "$_{0.30}$")
        name = name.replace("0.20", "$_{0.20}$").replace("0.10", "$_{0.10}$")
        name = name.replace("O3", "O$_3$")
    else:
        name = name.replace("2", "$_2$").replace("3", "$_3$").replace("4", "$_4$")
        name = name.replace("5", "$_5$").replace("x", "$_x$")
    return name

from uncertainties import ufloat

species = {
    "CO": {"C": 1, "O": 1},
    "CO2": {"C": 1, "O": 2},
    "CH4": {"C": 1, "H": 4},
    "C2H6": {"C": 2, "H": 6},
    "C2H4": {"C": 2, "H": 4},
    "C2H2": {"C": 2, "H": 2},
    "C3H8": {"C": 3, "H": 8},
    "C3H6": {"C": 3, "H": 6},
    "C3H4": {"C": 3, "H": 4},
    "C4H10": {"C": 4, "H": 10},
    "C3H4O2": {"C": 3, "H": 4, "O": 2},
    "CH3COOH": {"C": 2, "H": 4, "O": 2},
    "O2": {"O": 2},
    "H2": {"H": 2},
    "N2": {"N": 2},
    "H2O": {"H": 2, "O": 1}
}

species["methane"] = species["CH4"]
species["ethane"] = species["C2H6"]
species["ethene"] = species["C2H4"]
species["ethylene"] = species["C2H4"]
species["ethyne"] = species["C2H2"]
species["acetylene"] = species["C2H2"]
species["propane"] = species["C3H8"]
species["propene"] = species["C3H6"]
species["propylene"] = species["C3H6"]
species["butane"] = species["C4H10"]
species["acetic"] = species["C3H4O2"]
species["acrylic"] = species["CH3COOH"]
species["maleic"] = {"C": 4, "H": 2, "O": 3}

def p2XS(p, xC, xO, fuel="propane"):
    r = {}
    r["nCp"] = 0
    for s in p:
        if s != fuel:
            r["nCp"] += (ufloat(*p[s])/100) * species[s].get("C", 0)
    r["nCr"] = (ufloat(*p[fuel])/100) * species[fuel].get("C", 0)
    r["Xp"] = 100 * r["nCp"]/(xC/100 * species[fuel].get("C", 0))
    r["Xr"] = 100 * (xC/100 - r["nCr"]/species[fuel].get("C", 0))/(xC/100)
    r["Sp"] = {}
    for s in p:
        if s != fuel and species[s].get("C", 0) > 0:
            r["Sp"][s] = ufloat(*p[s]) * species[s].get("C", 0) / r["nCp"]
    r["nOr"] = (ufloat(*p["O2"])/100) * species["O2"].get("O", 0)
    r["XOr"] = 100 * (xO/100 - r["nOr"]/species["O2"].get("O", 0))/(xO/100)
    return r