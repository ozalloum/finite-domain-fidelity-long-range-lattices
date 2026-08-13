"""Core model and numerical utilities for finite-domain fidelity paper.

Model in rotating frame:
    i w_Z + L_alpha w + rho (|w|^2 - 1) w = 0.
Kernel:
    J_r = r^(-(1+alpha))/zeta(1+alpha).
"""
from __future__ import annotations
import numpy as np
import mpmath as mp
from scipy.special import zeta as hzeta
from scipy.optimize import brentq
from scipy.fft import next_fast_len

mp.mp.dps = 45


def zeta_norm(alpha: float) -> float:
    return float(hzeta(1.0 + alpha, 1.0))


def kernel(alpha: float, rmax: int) -> np.ndarray:
    r = np.arange(1, rmax + 1, dtype=float)
    return r ** (-(1.0 + alpha)) / zeta_norm(alpha)


def tail(alpha: float, m):
    """One-sided normalized tail sum sum_{r>m} J_r. m may be array-like."""
    return hzeta(1.0 + alpha, np.asarray(m) + 1.0) / zeta_norm(alpha)


def lambda_exact(alpha: float, q: float) -> float:
    """Infinite-lattice symbol Lambda(q) using polylogarithm."""
    s = 1 + mp.mpf(str(alpha))
    qq = mp.mpf(str(abs(float(q))))
    if qq == 0:
        return 0.0
    val = 2 * (1 - mp.re(mp.polylog(s, mp.e ** (1j * qq))) / mp.zeta(s))
    return float(val)


def lambda_symbol_float(alpha: float, q: float) -> float:
    """Double-precision evaluation of the same infinite polylogarithmic symbol.

    This is used for dense publication plotting grids; validation and selected-mode
    calculations continue to use ``lambda_exact``.
    """
    qq = abs(float(q))
    if qq == 0.0:
        return 0.0
    s = 1.0 + float(alpha)
    return float(2.0 * (1.0 - mp.fp.re(mp.fp.polylog(s, mp.fp.exp(1j * qq))) / mp.fp.zeta(s)))


def mi_gain(alpha: float, rho: float, q: float) -> float:
    lam = lambda_exact(alpha, q)
    return float(np.sqrt(max(lam * (2.0 * rho - lam), 0.0)))


def select_commensurate_mode(alpha: float, rho: float, periods=(64, 128, 256), min_ratio=0.95):
    """Find q* from Lambda(q*)=rho when possible and choose a nearby commensurate q_c."""
    f = lambda q: lambda_exact(alpha, q) - rho
    if f(np.pi) >= 0:
        qstar = brentq(f, 1e-10, np.pi, xtol=1e-12)
        gmax = rho
    else:
        qs = np.linspace(0, np.pi, 2001)
        gs = np.array([mi_gain(alpha, rho, q) for q in qs])
        i = int(np.argmax(gs))
        qstar, gmax = float(qs[i]), float(gs[i])
    choice = None
    for P in periods:
        m = max(1, int(round(qstar * P / (2 * np.pi))))
        qc = 2 * np.pi * m / P
        gc = mi_gain(alpha, rho, qc)
        ratio = gc / gmax if gmax > 0 else 0.0
        choice = (qstar, gmax, qc, gc, P, m, ratio, lambda_exact(alpha, qc))
        if ratio >= min_ratio:
            break
    return choice


class OpenOperator:
    """FFT-accelerated finite open operator and background-tail correction."""
    def __init__(self, alpha: float, N: int):
        self.alpha = float(alpha)
        self.N = int(N)
        J = kernel(alpha, N - 1)
        ker = np.r_[J[::-1], 0.0, J]
        self.nfft = next_fast_len(3 * N - 2)
        self.Kfft = np.fft.fft(ker, self.nfft)
        self.row_sum = self._conv(np.ones(N))
        idx = np.arange(N)
        self.kout = tail(alpha, idx) + tail(alpha, N - 1 - idx)

    def _conv(self, x):
        X = np.fft.fft(x, n=self.nfft, axis=-1)
        y = np.fft.ifft(X * self.Kfft, n=self.nfft, axis=-1)
        return y[..., self.N - 1:self.N - 1 + self.N]

    def linear_open(self, w):
        return self._conv(w) - self.row_sum * w

    def linear_corrected(self, w):
        # Background-tail closure: omitted exterior is held at w=1.
        return self.linear_open(w) + self.kout * (1.0 - w)

    def rhs_open(self, w, rho):
        return 1j * (self.linear_open(w) + rho * (np.abs(w) ** 2 - 1.0) * w)

    def rhs_corrected(self, w, rho):
        return 1j * (self.linear_corrected(w) + rho * (np.abs(w) ** 2 - 1.0) * w)


class PeriodicExactOperator:
    """Exact-periodized infinite operator at allowed Fourier modes.

    The positive-distance power-law weights are grouped by residue modulo N.
    Hurwitz zeta evaluates each residue class exactly, avoiding a large real-space cutoff.
    """
    def __init__(self, alpha: float, N: int):
        self.alpha = float(alpha)
        self.N = int(N)
        s = 1.0 + alpha
        zn = zeta_norm(alpha)
        W = np.zeros(N, dtype=float)
        # Residue 0 contains r=N,2N,... and does not affect Lambda because cos(0)=1.
        W[0] = N ** (-s)
        k = np.arange(1, N, dtype=float)
        W[1:] = (N ** (-s)) * hzeta(s, k / N) / zn
        # sum W = 1 (all positive r grouped by residue).
        self.lam = 2.0 * (1.0 - np.real(np.fft.fft(W)))
        self.lam[np.abs(self.lam) < 1e-14] = 0.0

    def linear(self, w):
        return np.fft.ifft(-self.lam * np.fft.fft(w, axis=-1), axis=-1)

    def rhs(self, w, rho):
        return 1j * (self.linear(w) + rho * (np.abs(w) ** 2 - 1.0) * w)


def rk4_step(w, h, rhs, rho):
    k1 = rhs(w, rho)
    k2 = rhs(w + 0.5 * h * k1, rho)
    k3 = rhs(w + 0.5 * h * k2, rho)
    k4 = rhs(w + h * k3, rho)
    return w + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def central_slice(N: int, fraction=0.25):
    width = max(4, int(round(N * fraction)))
    start = (N - width) // 2
    return slice(start, start + width)


def static_defects(alpha: float, N: int, q: float, central_fraction=0.25):
    n = np.arange(N)
    v = np.cos(q * n)
    lam = lambda_exact(alpha, q)
    Linf = -lam * v
    op = OpenOperator(alpha, N)
    Lo = np.real_if_close(op.linear_open(v))
    Lc = np.real_if_close(op.linear_open(v) - op.kout * v)  # perturbation form
    sl = central_slice(N, central_fraction)
    den = np.sqrt(np.mean(np.abs(v[sl]) ** 2))
    eo = np.sqrt(np.mean(np.abs(Lo[sl] - Linf[sl]) ** 2)) / den
    ec = np.sqrt(np.mean(np.abs(Lc[sl] - Linf[sl]) ** 2)) / den
    # Naive minimum-distance periodic closure.
    J = kernel(alpha, N // 2)
    r = np.arange(1, N // 2)
    lam_n = 2 * np.sum(J[:N // 2 - 1] * (1 - np.cos(r * q)))
    if N % 2 == 0:
        lam_n += J[N // 2 - 1] * (1 - np.cos((N // 2) * q))
    en = abs(lam_n - lam)
    return float(eo), float(ec), float(en)


def first_sustained_crossing(t, e, threshold, sustain=2):
    above = np.asarray(e) >= threshold
    if sustain <= 1:
        idx = np.where(above)[0]
        return float(t[idx[0]]) if len(idx) else np.nan
    for i in range(0, len(above) - sustain + 1):
        if np.all(above[i:i + sustain]):
            return float(t[i])
    return np.nan


def periodic_repeat(w_period, N):
    reps = int(np.ceil(N / len(w_period)))
    return np.tile(w_period, reps)[:N]
