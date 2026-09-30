"""
Worked examples of the error-mitigation section of proposal/quairk/ancillaStreaming/main.tex.

Every number quoted in the illustrative examples of that section is printed here.  Inputs taken from
the hardware tables of the paper are marked; everything else is a toy value chosen for illustration.

Usage:
    PYTHONPATH=. python3 -m qlb_port.mitigation_examples
"""

import json
import math


def main():
    print("1. shot noise of v = P(1) - P(0):  sigma = sqrt((1 - v^2) / N)")
    for N in (2000, 125):
        print(f"   N = {N:4d}, v = 0: sigma = {math.sqrt(1 / N):.4f}")

    print("\n2. coherent vs stochastic: |+> after N gates, each followed by RZ(eps); <X>")
    eps, N = 0.05, 16
    coh = math.cos(N * eps)
    p = math.sin(eps / 2) ** 2
    sto = (1 - 2 * p) ** N
    print(f"   eps = {eps}, N = {N}: coherent <X> = cos(N eps) = {coh:.3f} (error {1 - coh:.3f})")
    print(f"   twirled: p = sin^2(eps/2) = {p:.5f}, <X> = (1 - 2p)^N = {sto:.3f} (error {1 - sto:.3f})")

    print("\n3. single-qubit readout inversion, q1 of Table qcs-vel (e01 = 0.6%, e10 = 1.8%)")
    e01, e10, vm = 0.006, 0.018, 0.816
    print(f"   v_meas = {vm:+.3f} -> v = (v_meas - e01 + e10) / (1 - e01 - e10) = "
          f"{(vm - e01 + e10) / (1 - e01 - e10):+.3f}")

    print("\n4. post-selection toy: detected errors p_d, undetected p_u, erroneous shots average 0")
    pd, pu = 0.10, 0.05
    print(f"   p_d = {pd}, p_u = {pu}: raw = {1 - pd - pu:.3f} v, post-selected = "
          f"{(1 - pd - pu) / (1 - pd):.3f} v, shots kept {1 - pd:.2f}, "
          f"shot noise x {1 / math.sqrt(1 - pd):.3f}")

    print("\n5. phase sensitivity of v(t) = sin(2Et) to a phase shift phi = 0.1")
    E = float(json.load(open("qlb_port/qcs_quil_n2_cdr/qcs_results_cdr.json"))["E"])
    print(f"   2E = {2 * E:.4f}")
    for t in (1, 2, 3, 4):
        th = 2 * E * t
        print(f"   t = {t}: exact {math.sin(th):+.3f}, first order cos(2Et) phi = {math.cos(th) * 0.1:+.3f},"
              f" shift sin(2Et + phi) - sin(2Et) = {math.sin(th + 0.1) - math.sin(th):+.3f}")

    print("\n6. zero-noise extrapolation and rescaling with the values of Table qcs-zne")
    v1, v3 = 0.852, 0.623
    print(f"   t = 1 linear: v1 - (v3 - v1)/2 = {v1 - (v3 - v1) / 2:+.3f}")
    kap = math.log(v1 / v3) / 2
    print(f"   t = 1 two-point exponential: kappa = ln(v1/v3)/2 = {kap:.3f}, v0 = v1 e^kappa = {v1 * math.exp(kap):+.3f}")
    for t, v, m in ((1, 0.852, 0.550), (2, 0.178, 0.354), (3, -0.506, 0.293)):
        print(f"   t = {t} rescaled: v1 / sqrt(m) = {v:+.3f} / {math.sqrt(m):.3f} = {v / math.sqrt(m):+.3f}")

    print("\n7. Clifford data regression with the fit of Table qcs-cdr, t = 3")
    a, b, x = 0.031, 1.705, -0.550
    print(f"   a + b x = {a:+.3f} + {b:.3f} * ({x:+.3f}) = {a + b * x:+.3f}")
    print(f"   toy noise v_meas = 0.02 + 0.6 v_ideal inverts to b = 1/0.6 = {1 / 0.6:.3f}, a = -0.02/0.6 = {-0.02 / 0.6:+.4f}")

    print("\n8. repetition code with independent errors: P_L = sum_{k > d/2} C(d,k) p^k (1-p)^(d-k)")
    for p in (0.03, 0.2, 0.4):
        pl = [sum(math.comb(d, k) * p ** k * (1 - p) ** (d - k) for k in range(d // 2 + 1, d + 1))
              for d in (3, 5, 7, 9)]
        print(f"   p = {p}: d = 3, 5, 7, 9 -> " + ", ".join(f"{x:.2e}" for x in pl))


if __name__ == "__main__":
    main()
