"""
Physics validation for the corrected Quantum Lattice Boltzmann (QLB) Dirac solver.
==================================================================================
These tests validate the QLB scheme against INDEPENDENT exact references, not
against the solver's own output.  They are designed to catch the class of bug
that the previous test suite missed: an internally self-consistent scheme that
nonetheless does not solve the Dirac equation for a general spinor.

Representation (Dellar 2011): the three spatial streaming matrices are
{ALPHA_X, BETA, ALPHA_Z} and the mass/collision term acts along ALPHA_Y.

References used here:
  * Exact free-Dirac advection (massless): each ALPHA_X eigenmode translates at
    +/- c with no dispersion, so the QLB x-sweep must reproduce it to machine
    precision under the light-cone condition DX = c*DT.
  * Relativistic lattice dispersion: a plane wave e^{ikx} evolves with angular
    frequency |omega| = sqrt(sin^2 k + m^2); the QLB one-step operator must match
    this to second order in k.
  * Exact 2D Dellar-representation Dirac propagator
    H(kx,ky) = ALPHA_X sin(kx) + BETA sin(ky) + ALPHA_Y m,
    to which the operator-split QLB must converge as k -> 0.

Run:  python test_qlb_validation.py
"""

import numpy as np
from scipy.linalg import expm

import dirac_qlb_solver as s

ALPHA_X, ALPHA_Y, ALPHA_Z, BETA = s.ALPHA_X, s.ALPHA_Y, s.ALPHA_Z, s.BETA
CFG = s.SWEEP_CONFIG

_failures = []


def check(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        _failures.append(name)
    return ok


# ---------------------------------------------------------------------------
def collide(psi, cfg, mtil, gtil=0.0):
    """Characteristic-frame collision Q_char = a I - i b coll_gen (unitary)."""
    Om = mtil**2 - gtil**2
    D = 1 + Om / 4 - 1j * gtil
    a = (1 - Om / 4) / D
    b = mtil / D
    return a * psi - 1j * b * (psi @ cfg['coll_gen'].T)


def qlb_sweep_periodic(psi, axis, mtil, gtil=0.0):
    """One QLB sub-step on a periodic line/field using SWEEP_CONFIG (reference path)."""
    cfg = CFG[axis]
    dim = {'x': 0, 'y': 1, 'z': 2}[axis]
    pr = psi @ cfg['R_inv'].T
    pr = collide(pr, cfg, mtil, gtil)
    out = np.zeros_like(pr)
    for c in range(4):
        out[..., c] = np.roll(pr[..., c], cfg['stream'][c], axis=dim)
    return out @ cfg['R'].T


# ---------------------------------------------------------------------------
def test_matrix_properties():
    print("Test 1: rotation & collision matrix properties")
    I4 = np.eye(4)
    check("X unitary", np.allclose(s.X_INV_MATRIX @ s.X_MATRIX, I4))
    check("Z unitary", np.allclose(s.Z_INV_MATRIX @ s.Z_MATRIX, I4))
    check("X^-1 ALPHA_X X = diag(-1,-1,1,1)",
          np.allclose(s.X_INV_MATRIX @ ALPHA_X @ s.X_MATRIX, np.diag([-1., -1., 1., 1.])))
    check("Z^-1 ALPHA_Z Z = diag(-1,-1,1,1)",
          np.allclose(s.Z_INV_MATRIX @ ALPHA_Z @ s.Z_MATRIX, np.diag([-1., -1., 1., 1.])))
    check("streaming signs x=(-1,-1,1,1)", list(CFG['x']['stream']) == [-1, -1, 1, 1])
    check("streaming signs y=(1,1,-1,-1)", list(CFG['y']['stream']) == [1, 1, -1, -1])
    check("streaming signs z=(-1,-1,1,1)", list(CFG['z']['stream']) == [-1, -1, 1, 1])
    # {ALPHA_X, BETA, ALPHA_Z} mutually anticommute; mass ALPHA_Y anticommutes all
    def anti(A, B): return np.allclose(A @ B + B @ A, 0)
    check("{ALPHA_X,BETA,ALPHA_Z} mutually anticommute",
          anti(ALPHA_X, BETA) and anti(ALPHA_X, ALPHA_Z) and anti(BETA, ALPHA_Z))
    check("mass ALPHA_Y anticommutes streaming set",
          anti(ALPHA_Y, ALPHA_X) and anti(ALPHA_Y, BETA) and anti(ALPHA_Y, ALPHA_Z))
    # collision unitary with physical coefficients
    g, m = 0.3, 0.5
    Om = m*m - g*g; D = 1 + Om/4 - 1j*g; a = (1-Om/4)/D; b = m/D
    Qhat = a*np.eye(4) - 1j*b*ALPHA_Y
    check("Qhat = aI - i b ALPHA_Y unitary", np.allclose(Qhat.conj().T @ Qhat, np.eye(4)))


def test_massless_advection():
    print("Test 2: massless free x-sweep vs exact advection (generic spinor)")
    N = 256; x = np.arange(N)
    env = np.exp(-(x-80.)**2 / (2*10**2)) * np.exp(1j*0.3*x)
    sp = np.array([1, 0.4j, -0.3, 0.8]); sp = sp / np.linalg.norm(sp)
    psi0 = np.outer(env, sp)
    w, V = np.linalg.eigh(ALPHA_X)

    def exact(psi, t):
        coeff = psi @ V.conj()
        out = np.zeros_like(psi)
        for i in range(4):
            out += np.outer(np.roll(coeff[:, i], int(round(w[i]))*t), V[:, i])
        return out

    pq = psi0.copy()
    for _ in range(60):
        pq = qlb_sweep_periodic(pq, 'x', 0.0)
    pe = exact(psi0, 60)
    err = np.linalg.norm(pq - pe) / np.linalg.norm(pe)
    check("1D massless generic spinor matches exact advection", err < 1e-12, f"rel err {err:.1e}")
    check("norm conserved", abs(np.linalg.norm(pq)/np.linalg.norm(psi0) - 1) < 1e-12)


def test_dispersion():
    print("Test 3: massive dispersion vs sqrt(sin^2 k + m^2)")
    N = 256; x = np.arange(N)
    worst = 0.0
    for m in (0.1, 0.3):
        for k in (0.05, 0.1, 0.2):
            M = np.zeros((4, 4), complex)
            for c in range(4):
                pw = (np.exp(1j*k*x)[:, None] * np.eye(4)[c][None, :])
                outp = qlb_sweep_periodic(pw, 'x', m)
                M[:, c] = np.mean(outp / np.exp(1j*k*x)[:, None], axis=0)
            omega = np.sort(np.abs(np.angle(np.linalg.eigvals(M))))[-1]
            exact = np.sqrt(np.sin(k)**2 + m**2)
            worst = max(worst, abs(omega - exact))
    check("dispersion within 2nd-order tolerance", worst < 6e-3, f"max |domega| {worst:.1e}")


def test_2d_convergence():
    print("Test 4: 2D operator-split QLB converges to exact Dellar-rep Dirac")
    sp = np.array([1, 0.4j, -0.3, 0.8]); sp = sp / np.linalg.norm(sp)

    def exact2d(psi, m, T, N):
        pk = np.fft.fft2(psi, axes=(0, 1)); k = 2*np.pi*np.fft.fftfreq(N)
        out = np.zeros_like(pk)
        for a_ in range(N):
            for b_ in range(N):
                Hk = ALPHA_X*np.sin(k[a_]) + BETA*np.sin(k[b_]) + ALPHA_Y*m
                out[a_, b_] = expm(-1j*Hk*T) @ pk[a_, b_]
        return np.fft.ifft2(out, axes=(0, 1))

    errs = []
    for scale in (1, 2):
        N = 64*scale; x = np.arange(N); c = N//2; wid = 8*scale; k0 = 0.15/scale; T = 12; m = 0.1
        env = (np.exp(-((x[:, None]-c)**2 + (x[None, :]-c)**2)/(2*wid**2))
               * np.exp(1j*k0*(x[:, None] + 0.7*x[None, :])))
        psi = env[:, :, None] * sp[None, None, :]
        psi = psi[:, :, None, :]  # add NZ=1 axis for sweep API
        p0 = psi.copy()
        for _ in range(T):
            psi = qlb_sweep_periodic(psi, 'x', m/2)
            psi = qlb_sweep_periodic(psi, 'y', m/2)
        pe = exact2d(p0[:, :, 0, :], m, T, N)
        err = np.linalg.norm(psi[:, :, 0, :] - pe) / np.linalg.norm(pe)
        nrm = np.linalg.norm(psi) / np.linalg.norm(p0)
        errs.append(err)
        check(f"2D k0={k0:.3f}: err small & norm conserved",
              err < 0.12 and abs(nrm-1) < 1e-10, f"err {err:.2e}, norm {nrm:.5f}")
    check("2D error decreases as k halves", errs[1] < errs[0], f"{errs[0]:.2e} -> {errs[1]:.2e}")


def test_production_step_norm():
    print("Test 5: production perform_qlb_sub_step conserves probability (interior)")
    # massless free line, periodic-like interior check via norm on a padded line
    N = 200
    x = np.arange(N)
    env = np.exp(-(x-100.)**2/(2*12**2)) * np.exp(1j*0.25*x)
    char_right = np.array([0, 0, 1, 1], dtype=complex)/np.sqrt(2)
    sp = s.X_MATRIX @ char_right
    line = np.outer(env, sp)
    V = np.zeros(N)
    n0 = np.sum(np.abs(line)**2)
    out = line.copy()
    for _ in range(30):  # stays well inside the domain -> open BC drops nothing
        out = s.perform_qlb_sub_step(out, V, 'x', s.MATRICES, s.HBAR, s.C, 0.0, s.DT, s.DX)
    n1 = np.sum(np.abs(out)**2)
    # centre of mass should move +30 cells (pure right-mover)
    d0 = np.sum(np.abs(line)**2, axis=1); d1 = np.sum(np.abs(out)**2, axis=1)
    com0 = np.sum(x*d0)/np.sum(d0); com1 = np.sum(x*d1)/np.sum(d1)
    check("interior norm conserved", abs(n1/n0 - 1) < 1e-10, f"ratio {n1/n0:.8f}")
    check("physical +x propagation ~ +30 cells", abs((com1-com0) - 30) < 0.5, f"COM {com0:.1f}->{com1:.1f}")


if __name__ == "__main__":
    print("=" * 70)
    print("QLB Dirac solver — physics validation against exact references")
    print("=" * 70)
    test_matrix_properties()
    test_massless_advection()
    test_dispersion()
    test_2d_convergence()
    test_production_step_norm()
    print("=" * 70)
    if _failures:
        print(f"FAILED ({len(_failures)}): " + ", ".join(_failures))
        raise SystemExit(1)
    print("All physics validation tests passed.")
