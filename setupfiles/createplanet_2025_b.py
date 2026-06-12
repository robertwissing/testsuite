# PLEASE ERASE LINE ABOVE TO RUN THROUGH PYTHON3 INTERPRETER
# ERASE ---> %auto-ignore
# THE LINE WAS ADDED SO THE arXiV WOULD NOT PROCESS IT
###########################################################

# Python 3 code arXiV Version

# It needa SciPy and NumPy to run.
######################################################################
# Title: VARIABLE POLYPY, #CVS 1.06

# Author: Stephen Weppner, Eckerd College, March 2015. Version 1.0
# weppnesp@eckerd.edu

# This program calculates the interior of planets as detailed in
# A Variable Polytrope Index Applied to Planet and Material Models

# by S. P. Weppner, J. P. McKelvey, K. D. Thielen and A. K. Zielinski

# Submitted to the Monthly Notices of the Royal Astronomical Society
#######################################################################
# This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

########################################################
# > python3 supplemental_earth.py
#
# last lines of output:
#
# final radius, rho, pressure for layer  6
# 6382171.875 1800.00076139 67858.9821777
#
#
#   final mass       radius     density       gravity       pressure
# 5.9805557432e+24 6382171.875 5492.21393924 9.79896805169 364099932141.0
#
#  The full profile can be found in  earth_multiple.dat
#########################################################################

import sys
import scipy.integrate as integrate
import scipy.constants as constants
import numpy as np
import scipy.optimize as optimize
import mpmath

from planets_2025 import planet
from sklearn import svm
from varpoly import varpoly

# keep a few module-level constants/variables used by functions
pi = np.pi
G = constants.G
r_i00 = 500.00
hhmin = 1.e-9

# Forward-declared globals that legacy functions expect to exist at module level.
# These will be assigned inside run_model() before running the core solver.
material_layer = None
y1_0 = None
K_0 = None

#######################################################
# Helper to collect inputs (wrap original input section)
#######################################################

def get_inputs(planet_name="testplanet_2layer"):
    """Return a dict with all input variables produced by the original input section.

    Call with get_inputs("testplanet_2layer") or another planet name as in planets_2025.py
    """
    pl = planet(planet_name)
    file_string = pl[0]
    p_00 = pl[1]
    Tc = pl[2]
    r_end = pl[3]
    Tcond = pl[4]
    materials = pl[5]
    num_of_materials = len(materials)

    matindex = list(range(num_of_materials))
    yy1_0 = [m.rho0 for m in materials]
    BB_0 = [m.B0 for m in materials]
    nn = [m.dBdP0 for m in materials]
    ZZ = [m.ZZ for m in materials]
    AA = [m.AA for m in materials]
    tminl = [m.tmin for m in materials]
    deltal = [m.delta for m in materials]
    a_exp = [m.Pexp_const_a for m in materials]
    b_exp = [m.Pexp_const_b for m in materials]
    c_exp = [m.grun_const_c for m in materials]

    # Calculate Cv (Dulong-Petit limit) for all layers - Cv=3R/avg(M)
    Cvinf = np.zeros(num_of_materials)
    for j in range(num_of_materials):
        Cvinf[j] = 24943 / (AA[j])

    return dict(
        file_string=file_string,
        p_00=p_00,
        Tc=Tc,
        r_end=r_end,
        Tcond=Tcond,
        materials=materials,
        num_of_materials=num_of_materials,
        matindex=matindex,
        yy1_0=yy1_0,
        BB_0=BB_0,
        nn=nn,
        ZZ=ZZ,
        AA=AA,
        tminl=tminl,
        deltal=deltal,
        a_exp=a_exp,
        b_exp=b_exp,
        c_exp=c_exp,
        Cvinf=Cvinf,
    )


##########################################
# Core functions (kept largely unchanged)
##########################################

def deriv_first(u, r):
    # u1 = r * rho, u2 = du1/dr, u3 = m, u4 = pressure, ut = temperature
    u1, u2, u3, u4, ut = u
    u1_deriv = u2  # du/dr
    # u = r * rho for small r version
    [n_dyn, B, P, comp_const, gamma] = calculate_variables(u1 / r, ut)

    # d^2 u/dr^2
    u2_deriv = -4.0 * pi * G * \
        (comp_const * y1_0**n_dyn * K_0) * u1**(3 - n_dyn) * r**(n_dyn - 2)
    u3_deriv = 4 * pi * r * u1  # dm/dr
    u4_deriv = u1 * u3 * G / r**3    # dp/dr
    ut_deriv = -((gamma * ut) / B) * u4_deriv

    return [u1_deriv, u2_deriv, u3_deriv, u4_deriv, ut_deriv]


def deriv_y(y, r):
    # y1 = rho, y2 = d rho/dr, y3 = m, y4 = pressure
    y1, y2, y3, y4, yt = y

    [n_dyn, B, P, comp_const, gamma] = calculate_variables(y1, yt)
    y1_deriv = y2  # drho/dr
    y2_deriv = -((n_dyn - 2) / y1 * y2**2 + 2.0 / r * y2 + 4.0 * pi * G * comp_const
                 * y1_0**n_dyn * K_0 * y1**(3 - n_dyn))  # d^2rho/dr^2
    y3_deriv = 4 * pi * r**2 * y1   # dm/dr
    y4_deriv = y1 * y3 * G / r**2    # dp/dr
    yt_deriv = -((gamma * yt) / B) * y4_deriv
    if abs(n_dyn - calculate_n(y1, yt, n_dyn, B)) > .00001:
        print("N Comparison failure!", calculate_n(y1, yt, n_dyn, B), n_dyn)
    comp_const2 = (y1**n_dyn / B) / (y1_0**n_dyn * K_0)
    if abs(comp_const - comp_const2) > .00001:
        print("Constant Comparison failure! ", comp_const, comp_const2, " r ", r, " density ", y1, " rho0 ", material_layer.rho0, " likely in expanded state, reduce central pressure")
    return [y1_deriv, y2_deriv, y3_deriv, y4_deriv, yt_deriv]


def calculate_variables(y1, yt):
    if(y1 <= material_layer.rho0):
        rho = material_layer.rho0
    else:
        rho = y1
    material_layer.calculate_P(y1,yt)
    material_layer.calculate_B(y1,yt)
    material_layer.calculate_dBdP(y1,yt)
    material_layer.calculate_grun_cold(y1)
    B = material_layer.B[0]
    nad = material_layer.dBdP[0]
    P = material_layer.P[0]
    gamma = material_layer.grun[0]
    comp_const = calculate_comp_const(y1, yt, nad, B)  # weighting function
    return [nad, B, P, comp_const, gamma]


def integrand(b, n, x):
    return mpmath.exp(-x * b) / b**n


def integrandU(b, T):
    Pt0 = calculate_PisoT0(b)
    return Pt0 / b**2


def expn(n, x):
    return integrate.quad(integrand, 1, np.inf, args=(n, x))[0]


def func_rho_P(rho, p_0, T):
    [n_dyn, B, P, comp_const, gamma] = calculate_variables(rho, T)
    P = np.array(P, dtype=float)
    return P - p_0


def func_rho_P_deriv(rho, p_0, T):
    [n_dyn, B, P, comp_const, gamma] = calculate_variables(rho, T)
    B = np.array(B, dtype=float)
    rho = np.array(rho, dtype=float)
    return B / rho


def calculate_niso(rho):
    if(rho <= y1_0):
        rho = y1_0
    B_core = calculate_B_coreisoT0(rho)
    K_core = 1.0 / B_core
    lncomp_const = mpmath.log(calculate_comp_constisoT0(rho))
    n_core = (lncomp_const + mpmath.log(K_0 / K_core)) / mpmath.log(rho / y1_0)
    return n_core


def calculate_n(rho, T, n, B_core):
    K_core = 1.0 / B_core
    lncomp_const = mpmath.log(calculate_comp_const(rho, T, n, B_core))
    n_core = (lncomp_const + mpmath.log(K_0 / K_core)) / mpmath.log(rho / y1_0)
    return n_core


def calculate_comp_constisoT0(rho):
    if(rho <= y1_0):
        rho = y1_0
    B_0 = 1 / K_0
    lncomp_const = ((calculate_dBdpisoT0(rho)) * mpmath.log(rho / y1_0)
                    + mpmath.log(B_0 / calculate_B_coreisoT0(rho)))
    return np.exp(lncomp_const)


def calculate_comp_const(rho, T, n, B):
    lncomp_const = (n * mpmath.log(rho / material_layer.rho0)
                    + mpmath.log(material_layer.B0 / B))
    return mpmath.exp(lncomp_const)


# This is the main function, it is called for each layer going from low r to high r.
# It solves the Lane-Emden Equation as defined in the paper.
def calculate_compound_planet(rho_0, r_start, g_start, m_start, p_start, T_start, material_index):
    global complete, num_final
    # initialize variables
    # ensure functions see the intended module-level variables; the main caller sets them
    material_layer = globals().get('material_layer')
    complete = False
    num_final = 0
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    yt = []
    r = []
    # conditions to change step size which is defined by a rho/rho_0 argument.
    r_i0 = r_i00
    if rho_0 < 1.5 * y1_0:
        r_i0 *= 0.5
    if rho_0 < 1.25 * y1_0:
        r_i0 *= 0.5
    if rho_0 < 1.10 * y1_0:
        r_i0 *= 0.5
    if rho_0 < 1.05 * y1_0:
        r_i0 *= 0.25
    # now special case for core, r = 0 (first step only)
    if material_index == 0:
        # use approximation for first 11 steps
        # (in 11 steps, 0 to r_i0 meters)
        r_0 = np.linspace(1, r_i0, 11)
        # intial values for u = r*rho, u', m, P
        # the initial values assume r = 1 meter.
        uinit = [rho_0, rho_0, 0.0, 0.0, T_start]
        u = integrate.odeint(deriv_first, uinit, r_0, rtol=1e-9, mxhnil=0)
        u1, u2, u3, u4, ut = np.hsplit(u, 5)
        u1 = np.reshape(u1, 11)
        u2 = np.reshape(u2, 11)
        u3 = np.reshape(u3, 11)
        u4 = np.reshape(u4, 11)
        ut = np.reshape(ut, 11)
        r.extend([0.0])
        y1.extend([rho_0])
        y2.extend([0.0])
        y3.extend([0.0])
        y4.extend([0.0])
        yt.extend([T_start])
        # pressure not needed so actually calculating the
        # difference from the core pressure
        u1_f = u1[-1] / r_0[-1]   # in terms of y
        u2_f = u2[-1] / r_0[-1] - u1[-1] / r_0[-1] / r_0[-1]
        u3_f = u3[-1]
        u4_f = u4[-1]
        ut_f = ut[-1]
        r.extend([r_i0 * 1.0])
        y1.extend([u1_f])
        y2.extend([u2_f])
        y3.extend([u3_f])
        y4.extend([u4_f])
        yt.extend([ut_f])
        yinit = [u1_f, u2_f, u3_f, u4_f, ut_f]
        r_f = r_i0 + 20 * r_i0
        r_i = r_i0
    else:  # this section only sets up initial values in the 2nd layer and beyond
        # the true work is done in the while loop below
        r_i = r_start
        r_f = r_start + 20 * r_i0
        # found analytically
        [n_start, B_start, P_start2, comp_const_start,
            gamma_start] = calculate_variables(rho_0, T_start)
        # This is the derivative constraint found in the paper
        deriv_0 = (-y1_0**n_start * K_0 / (rho_0**(n_start - 2))
                   * g_start * comp_const_start)
        r.extend([r_start])
        y1.extend([rho_0])
        y2.extend([deriv_0])
        y3.extend([m_start])
        y4.extend([p_start])
        yt.extend([T_start])
        # pressure not needed any more so actually calculating
        # the relative difference from core pressure
        # the moment where we set the initial values for layer 2 and beyond
        yinit = [rho_0, deriv_0, m_start, p_start, T_start]
    reach_final = True
    # This loop below is where most of the work is done -- solving the D.E.Q.s
    while yinit[0] > y1_0 and reach_final:
        rr = np.linspace(r_i, r_f, 21)
        y, info = integrate.odeint(deriv_y, yinit, rr, rtol=1.e-9, full_output=1,
                                   hmin=hhmin, mxstep=5000000, mxhnil=0)
        temp1, temp2, temp3, temp4, tempT = np.hsplit(y, 5)
        temp1 = np.delete(temp1, 0)
        temp2 = np.delete(temp2, 0)
        temp3 = np.delete(temp3, 0)  # delete duplicate points
        temp4 = np.delete(temp4, 0)
        tempT = np.delete(tempT, 0)
        rr = np.delete(rr, 0)
        y1.extend(temp1)   # add to original dataset
        y2.extend(temp2)
        y3.extend(temp3)
        y4.extend(temp4)
        yt.extend(tempT)
        r.extend(rr)
        # these blocks search for surface of planet from full set, tested
        num = 20
        found = False
        l1 = len(temp1)
        while num > 0 and not found:
            if info['hu'][l1 - num] < hhmin:
                found = True
                reach_final = False
                if temp1[l1 - num] <= y1_0:
                    if num_of_materials - 1 != material_index:
                        print("ERROR!, interior material hit surface pressure for layer ",
                              material_index + 1, ", shrink radius")
                        sys.exit(1)
                    else:
                        complete = True
            # end of layer check
            if rr[l1 - num] > r_end[material_index]:
                found = True
                reach_final = False
                complete = True
            num -= 1
            num_final = num
        if reach_final:
            r_i = r_f
            r_f = r_i + 20 * r_i0
            yinit = [temp1[-1], temp2[-1], temp3[-1], temp4[-1], tempT[-1]]
    # delete excess points beyond layer boundary
    ii = np.arange(len(y1) - num_final - 1, len(y1), 1)
    y1 = np.delete(y1, ii)
    y2 = np.delete(y2, ii)
    y3 = np.delete(y3, ii)
    y4 = np.delete(y4, ii)
    yt = np.delete(yt, ii)
    r = np.delete(r, ii)
    y = np.array([r, y1, y2, y3, y4, yt])
    rho_last = y1[-1]
    [n_last, B_last, P_last, comp_const_last,
        gamma_last] = calculate_variables(y1[-1], yt[-1])
    g_last = y3[-1] * G / r[-1] / r[-1]
    # we can check the accuracy of D.E.Q. by comparing it to analytical
    # version of the derivative which appears in the paper
    d_test = (-y1_0**n_last * K_0 / (rho_last**(n_last - 2))
              * g_last * comp_const_last)
    print("DERIVATIVE CHECK AT BOUNDARY", d_test, y2[-1])
    return y


############################################################
# run_model() wraps the original top-level main program. The
# legacy functions above expect several names to exist in the
# module global namespace (material_layer, y1_0, K_0, etc.).
# run_model assigns those module globals before calling the
# solver routines so minimal changes to the original code are
# required.
############################################################

def run_model(inputs):
    """Run the planet solver using the inputs returned by get_inputs().

    Returns a dict with final results and arrays.
    """
    # declare module globals that the legacy functions expect
    global material_layer, y1_0, K_0, pi, G

    # unpack inputs
    file_string = inputs['file_string']
    p_00 = inputs['p_00']
    Tc = inputs['Tc']
    r_end = inputs['r_end']
    Tcond = inputs['Tcond']
    materials = inputs['materials']
    num_of_materials = inputs['num_of_materials']
    matindex = inputs['matindex']
    yy1_0 = inputs['yy1_0']
    BB_0 = inputs['BB_0']
    nn = inputs['nn']
    ZZ = inputs['ZZ']
    AA = inputs['AA']
    tminl = inputs['tminl']
    deltal = inputs['deltal']
    a_exp = inputs['a_exp']
    b_exp = inputs['b_exp']
    c_exp = inputs['c_exp']
    Cvinf = inputs['Cvinf']

    globals()['r_end'] = r_end
    globals()['num_of_materials'] = num_of_materials
    globals()['Tcond'] = Tcond
    globals()['p_00'] = p_00
    
    # initialize variables that used to be top-level
    r_start = 0.0
    g_start = 0.0
    m_start = 0.0
    p_start = 0.0
    T_start = Tc
    p_0 = p_00
    rho_start = 1.e9
    print("PolyPy -- version 1")
    materiallist = np.array([])

    # prepare arrays used in the top-level loop
    nn_rho = []
    nn_press = []
    nn_mass = []
    nn_radius = []
    p_theory = []
    super_y = np.empty(5)
    y5 = []
    y6 = []
    y7 = []
    y8 = []
    y9 = []
    y10 = []
    y11 = []
    y12 = []
    y13 = []
    y14 = []
    y15 = []
    y16 = []

    # main loop over layers (kept as in original)
    for index1 in range(num_of_materials):
        material_index = index1
        # set module-level material_layer so calculate_variables and other
        # legacy routines see the correct material.
        material_layer = materials[material_index]
        globals()['material_layer'] = material_layer

        print("\n")
        print("****** Layer ", material_index + 1, " ******")
        n0 = material_layer.dBdP0
        y1_0 = material_layer.rho0
        globals()['y1_0'] = y1_0
        K_0 = 1.0 / material_layer.B0
        globals()['K_0'] = K_0
        rho_critical = material_layer.rho_critical
        print("COMPARE B", material_layer.calculate_B_cold(rho_critical + .01),
              material_layer.calculate_B_cold(rho_critical - .01))
        print("COMPARE dBdp", rho_critical, material_layer.calculate_dBdP_cold(rho_critical + .001),
              material_layer.calculate_dBdP_cold(rho_critical - .001))
        chi_tot = 0.0
        chi_array = np.array([])
        A_0 = material_layer.A0
        A_1 = material_layer.A1
        A_2 = material_layer.A2
        print("           n   rho_0   K_0       A_0            A_1           A_2 ")
        print("Starting ", n0, y1_0, K_0, A_0, A_1, A_2)
        print("Starting P and T is ", p_0/10**9, T_start)
        print("1-A1+A2/A1 ",1-(A_1+A_2)/A_1)
        print("A0/A1", A_0/A_1)
        print("expint", expn(1.1, 0.01))
        rho_0 = optimize.newton(func_rho_P, y1_0*3.0,func_rho_P_deriv, args=(p_0, T_start,))
        print("Rho initial found:", rho_0, "  pressure:", p_0/10**9)
        if rho_0 > rho_start:
            print("Starting density for interior material is too low")

        print("Checking boundaries:")
        material_layer.calculate_P_cold(y1_0)
        material_layer.calculate_B_cold(y1_0)
        material_layer.calculate_dBdP_cold(y1_0)
        material_layer.calculate_grun_cold(y1_0)

        Bi0 = material_layer.B_c[0]
        ni0 = material_layer.dBdP_c[0]
        n2i0 = material_layer.d2BdP_c[0]
        n3i0 = material_layer.d3BdP_c[0]
        Pi0 = material_layer.P_c[0]
        [ti0, dti0, d2ti0, d3ti0] = [material_layer.t[0], material_layer.dtdP[0], material_layer.d2tdP[0] , material_layer.d3tdP[0] ]
        [gammai0, qi0, qlambi0] = [material_layer.grun[0],material_layer.grun_q[0],material_layer.grun_lamb[0]]

        rhoinf = 10**15
        material_layer.calculate_P_cold(rhoinf)
        material_layer.calculate_B_cold(rhoinf)
        material_layer.calculate_dBdP_cold(rhoinf)
        material_layer.calculate_grun_cold(rhoinf)
        Btinf = material_layer.B_c[0]
        ntinf = material_layer.dBdP_c[0]
        n2tinf = material_layer.d2BdP_c[0]
        n3tinf = material_layer.d3BdP_c[0]
        Ptinf = material_layer.P_c[0]
        [tinf, dtinf, d2tinf, d3tinf] = [material_layer.t[0], material_layer.dtdP[0], material_layer.d2tdP[0] , material_layer.d3tdP[0]]
        [gammainf, qinf, qlambinf] = [material_layer.grun[0],material_layer.grun_q[0],material_layer.grun_lamb[0]]

        material_layer.calculate_P_cold(rho_0)
        material_layer.calculate_B_cold(rho_0)
        material_layer.calculate_dBdP_cold(rho_0)
        material_layer.calculate_grun_cold(rho_0)
        Bis = material_layer.B_c[0]
        nis = material_layer.dBdP_c[0]
        n2is = material_layer.d2BdP_c[0]
        n3is = material_layer.d3BdP_c[0]
        Pis = material_layer.P_c[0]
        [tis, dts, d2ts, d3ts] = [material_layer.t[0], material_layer.dtdP[0], material_layer.d2tdP[0] , material_layer.d3tdP[0] ]
        [gammas, qs, qlambs] = [material_layer.grun[0],material_layer.grun_q[0],material_layer.grun_lamb[0]]

        print("Pi0: ", Pi0/10**9, "Pinf: ", Ptinf/10**9)
        print("Bi0: ", Bi0, "Binf: ", Btinf)
        print("ni0: ", ni0, "ninf: ", ntinf)
        print("n2i0: ", n2i0, "n2inf: ", n2tinf)
        print("n3i0: ", n3i0, "n3inf: ", n3tinf)
        print("ti0: ", ti0, "tinf: ", tinf)
        print("dti0: ", dti0, "dtinf: ", dtinf)
        print("d2ti0: ", d2ti0, "d2tinf: ", d2tinf)
        print("d3ti0: ", d3ti0, "d3tinf: ", d3tinf)
        print("gammai0: ", gammai0, "gammainf: ", gammainf)
        print("qi0: ", qi0, "qinf: ", qinf)
        print("qlambi0: ", qlambi0 / qi0, "lambinf: ", qlambinf / qinf)
        print("gamma: ", gammas)
        print("q: ", qs)
        print("t: ", tis)
        print("asymtotic relation ninf*Pinf/Binf = 1 ", ntinf*Ptinf/Btinf)
        print("asymtotic relation for lambinf ", (Btinf*n2tinf/(1-ntinf*Ptinf/Btinf))/ntinf+ntinf )
        print("relation between ninf and lambinf should be -(ntinf+lambinf) = -2 ",Btinf**2*n3tinf/(Btinf*n2tinf))

        # This calls the routine which solves the Lane - Emden Eq.
        y = calculate_compound_planet(
            rho_0, r_start, g_start, m_start, p_start, T_start, index1)
        if material_index == 0:
            super_y = y
        else:
            super_y = np.hstack((super_y, y))

        r, y1, y2, y3, y4, yt = np.vsplit(y, 6)
        si = r.shape[1]
        r = np.reshape(r, si)
        y1 = np.reshape(y1, si)
        y2 = np.reshape(y2, si)
        y3 = np.reshape(y3, si)
        y4 = np.reshape(y4, si)
        yt = np.reshape(yt, si)
        for index5 in range(len(y1)):
            myrho=y1[index5]
            mytemp=yt[index5]
            material_layer.calculate_P(myrho,mytemp)
            material_layer.calculate_B(myrho,mytemp)
            material_layer.calculate_dBdP(myrho,mytemp)
            material_layer.calculate_grun_cold(myrho)
            material_layer.calculate_U(myrho,mytemp)
            P = material_layer.P[0]
            Pt0=material_layer.P_c[0]
            Pth=material_layer.P_vib[0]
            B=material_layer.B[0]
            U=material_layer.U[0]
            n_dyn = material_layer.dBdP[0]
            n_iso = material_layer.dBdP_iso[0]
            gamma = material_layer.grun[0]
            q = material_layer.grun_q[0]
            qlamb = material_layer.grun_lamb[0]
            t = material_layer.t[0]
            Cv = material_layer.Cv[0]
            print(Pth/P);
            y5.extend([n_dyn])
            y7.extend([B])
            y8.extend([Pth])
            y9.extend([Pth/P])
            y10.extend([gamma])
            y11.extend([q])
            y12.extend([qlamb])
            y13.extend([t])
            y14.extend([Cv])
            y15.extend([U])
            y16.extend([P])
            if index5 == 0:
                y6.extend([0.0])
            else:
                y6.extend([G * y3[index5] / (r[index5])**2])

        # do gravity calculation
        voly = np.array(r) * np.array(r) * np.array(y1)
        gr = integrate.cumtrapz(y=voly, x=r, initial=0.0)
        nn_radius.extend([r[-1]])
        nn_mass.extend([gr[-1] * 4 * pi])
        gr = gr * 4 * pi * G
        gr[0] = 0.0
        for num in range(1, len(gr)):
            gr[num] = gr[num] / r[num] / r[num]
        # do pressure
        voly2 = np.array(y1) * np.array(gr)
        p = integrate.simps(y=voly2, x=r)
        p_theory = np.append(p_theory, mpmath.log10(p))
        r_start = r[-1]
        g_start = y3[-1] * G / (r[-1])**2
        m_start = y3[-1]
        p_start = y4[-1]
        rho_start = y1[-1]
        T_start = yt[-1]
        p_0 = p_00 - p_start  # change starting pressure
        if(Tcond[index1 + 1] > 0):
            T_start = Tcond[index1 + 1]
        print("final radius, rho, pressure, temperature and mass for layer(earth mass) ", material_index + 1)
        print(r[-1], y1[-1], p_0/10**9, yt[-1],y3[-1]/(5.972*10**24))
        materiallayer = np.zeros(r.size)
        materiallayer.fill(material_index)
        materiallist = np.append(materiallist, materiallayer)
        # End loop over layers

    # reshape for final print out
    r, y1, y2, y3, y4, yt = np.vsplit(super_y, 6)
    si = r.shape[1]
    r = np.reshape(r, si)
    y1 = np.reshape(y1, si)
    y2 = np.reshape(y2, si)
    y3 = np.reshape(y3, si)
    y4 = np.reshape(y4, si)
    yt = np.reshape(yt, si)
    data_out = np.column_stack((r, y1, y2, y3, y6, y16, y8, yt, materiallist, y5,
                                y10, y11, y12, y13, y14, y9,
                                y15))
    np.savetxt(file_string, data_out, fmt='%1.8e',
               header='[r(m)] [rho(kg/m^3)] [drho/dr(kg/m^4)] [mass(kg)] [g(m/s^2)] [P (Pa)] [Pth] [temperature] [material] [dB/dP] [grun] [grun_q] [grun_qlamb] [t-param] [C_v] [Pth/P] [U] ')
    volume = (4. / 3.) * pi * r[-1]**3

    data_out_material = np.column_stack(( matindex, yy1_0, BB_0, nn, AA, ZZ, tminl, deltal, a_exp, b_exp, c_exp))
    np.savetxt(file_string + "_material", data_out_material, fmt='%1.4e',
               header='[material] [rho_0(kg/m^3)] [B_0 ] [dB/dP_0] [AA] [ZZ] [t_min] [delta] [a_exp] [b_exp] [c_exp]')


    print("\n")
    print("  final mass in earth masses       radius in earth radius radius in meter    average density     gravity on surface   pressure on surface      temperature")
    print(y3[-1]/(5.972*10**24), r[-1]/(6371*10**3), r[-1], y3[-1] / volume, y6[-1], y16[-1]/10**9, yt[-1])

    print('\n', "The full profile can be found in ", file_string)

    results = dict(
        final_mass=y3[-1],
        final_radius=r[-1],
        final_density=y3[-1]/volume,
        final_gravity=y6[-1],
        final_pressure=y16[-1],
        final_temperature=yt[-1],
        volume=volume,
        r=r,
        rho=y1,
        materiallist=materiallist,
        file_string=file_string,
    )
    return r,y1,yt,materiallist,y3[-1],r[-1]  ##return r,dens,temp,final mass(kg), final rad(m)


if __name__ == "__main__":
    inputs = get_inputs("testplanet_2layer")
    results = run_model(inputs)
    print("\nFinal Results:")
    for k, v in results.items():
        if isinstance(v, (int, float)):
            print(f"{k:20s} = {v}")



import os
import numpy as np

def run_or_load(planet_name="testplanet_2layer", force=False):
    """
    If output file for `planet_name` exists and force is False, load it and
    return (r, rho, T, final_mass, final_radius) to match run_model()'s tuple.
    Otherwise call run_model(inputs) and return its result.
    """
    inputs = get_inputs(planet_name)
    file_string = inputs['file_string']

    # If file exists and not forcing a re-run -> load and return same tuple shape
    if (not force) and os.path.exists(file_string):
        print(f"Found existing profile '{file_string}' — loading instead of running.")
        data = np.loadtxt(file_string)   # np.savetxt wrote header as comments, loadtxt skips them
        if data.ndim == 1:
            # single-line file -> make it a 2D array with one row
            data = data.reshape(1, -1)

        # Column mapping based on your np.savetxt(..., column_stack(...)):
        # 0: r, 1: rho, 2: drho/dr, 3: mass, 4: g, 5: P, 6: dB/dP, 7: temperature, 8: material, ...
        try:
            r = data[:, 0]
            rho = data[:, 1]
            temp = data[:, 7]
            material = data[:, 8]
            mass_arr = data[:, 3]
        except IndexError:
            raise RuntimeError(f"Loaded file {file_string} doesn't have expected columns (need >=8 columns).")

        final_mass = mass_arr[-1]
        final_radius = r[-1]
        return r, rho, temp, material, final_mass, final_radius

    # Otherwise run the model
    print(f"No existing profile found for '{file_string}' or force=True -> running model.")
    res = run_model(inputs)

    # run_model might already return the expected tuple (r,rho,temp,final_mass,final_rad).
    # Accept that directly; if it returns a dict, try to extract the fields.
    if isinstance(res, tuple) or isinstance(res, list):
        return tuple(res)
    elif isinstance(res, dict):
        # try to map dictionary keys to the expected order
        r = res.get('r')
        rho = res.get('rho') or res.get('y1')
        temp = res.get('yt') or res.get('temperature') or res.get('T')
        material = res.get('materiallist') or res.get('material')
        final_mass = res.get('final_mass') or (res.get('y3')[-1] if res.get('y3') is not None else None)
        final_radius = res.get('final_radius') or (res.get('r')[-1] if res.get('r') is not None else None)

        if any(v is None for v in (r, rho, temp, material, final_mass, final_radius)):
            raise RuntimeError("run_model returned a dict but missing fields needed to form the tuple.")
        return r, rho, temp, material, final_mass, final_radius
    else:
        raise RuntimeError("run_model returned an unexpected type. Expected tuple or dict.")

