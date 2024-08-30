import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import CoolProp as CP

"""
Conversion between fractions:
https://www.researchgate.net/file.PostFileLoader.html?id=59422cefed99e1792a3437e0&assetKey=AS%3A505403965480960%401497509103277
"""

fontsize = 15
mpl.rcParams.update({"font.family": "Times New Roman",
                     "mathtext.fontset": "dejavuserif",
                     "font.size": fontsize,
                     "axes.linewidth": 1.,
                     "axes.labelsize": fontsize,
                     "axes.titlesize" : fontsize,
                     "legend.fontsize": fontsize,
                     "xtick.top": True,
                     "xtick.bottom": True,
                     "ytick.left": True,
                     "ytick.right": True,
                     "xtick.direction": "in",
                     "ytick.direction": "in",
                     "xtick.major.pad": 7,
                     "ytick.major.pad": 7})

HEOS = {"air": CP.AbstractState("HEOS", "Air"),
        "ch4": CP.AbstractState("HEOS", "CH4"),
        "co2": CP.AbstractState("HEOS", "CO2")}

def get_mole_vol_fractions(T_degC, p_bara):
    for k in HEOS.keys():
        HEOS[k].update(CP.PT_INPUTS, p_bara*1e5, T_degC+273.15)
    # mole fraction
    x_ch4 = (HEOS["air"].molar_mass() / HEOS["co2"].molar_mass() - 1) / (HEOS["ch4"].molar_mass() / HEOS["co2"].molar_mass() - 1)
    # vol fraction
    phi_ch4 = HEOS["ch4"].compressibility_factor() * x_ch4 / (HEOS["ch4"].compressibility_factor() * x_ch4 + HEOS["co2"].compressibility_factor() * (1 - x_ch4))
    return {"x_ch4": x_ch4, "phi_ch4": phi_ch4}

fig, ax = plt.subplots(1, 1, figsize=(5, 4))

T = np.arange(20, 25, 1)
p = np.arange(1, 2.1, 0.1)
for _T in T:
    phis = [get_mole_vol_fractions(_T, _p)['phi_ch4'] for _p in p]
    ax.plot(p, phis, label=str(_T)+r"$~\degree{\rm C}$")
ax.set_xlabel(r"$p~/~{\rm bar(a)}$")
ax.set_ylabel(r"$\phi_{\rm CH4}$")
ax.set_xlim((1, 2))
ax.legend(frameon=False)
fig.tight_layout(pad=0.2)
plt.show()