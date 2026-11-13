"""Constantes y fórmulas comunes. Todo en SI."""
import numpy as np

KB = 1.380649e-23          # J/K, exacto (SI 2019)
T_AMBIENTE = 296.15        # K (23 °C), valor por defecto si la grabación no trae temperatura


def gamma_stokes(eta, R):
    """Arrastre de Stokes de una esfera lejos de las paredes.

    eta: viscosidad [Pa·s]; R: radio [m]. Devuelve gamma = 6 pi eta R [N·s/m].
    """
    return 6.0 * np.pi * eta * R


def faxen_paralelo(R, h):
    """Factor de Faxén para el arrastre paralelo a una pared plana: gamma_par / gamma_stokes.

    R: radio [m]; h: distancia del CENTRO de la bead a la pared [m] (h > R).
    gamma_par = 6 pi eta R / (1 - 9/16 (R/h) + 1/8 (R/h)^3 - 45/256 (R/h)^4 - 1/16 (R/h)^5)
    """
    q = np.asarray(R, float) / np.asarray(h, float)
    return 1.0 / (1 - 9 / 16 * q + q ** 3 / 8 - 45 / 256 * q ** 4 - q ** 5 / 16)


def viscosidad_glicerina_agua(fraccion_masa, T_K):
    """Viscosidad de una mezcla glicerina-agua [Pa·s]. Cheng, Ind. Eng. Chem. Res. 47, 3285 (2008).

    fraccion_masa: masa de glicerina / masa total (0..1). T_K: temperatura [K] (vale 0-100 °C).
      mu_w = 1.790 exp((-1230 - T) T / (36100 + 360 T))   [mPa·s], T en °C
      mu_g = 12100 exp((-1233 + T) T / (9900 + 70 T))
      a = 0.705 - 0.0017 T,  b = (4.9 + 0.036 T) a^2.5
      alpha = 1 - C + a b C (1 - C) / (a C + b (1 - C)),   mu = mu_w^alpha mu_g^(1 - alpha)
    Da ~1.76 mPa·s con 20 % en masa a 20 °C. A VALIDAR por Lenny (¿el 20 % es en masa o en volumen?).
    """
    T = np.asarray(T_K, float) - 273.15
    C = np.asarray(fraccion_masa, float)
    mu_w = 1.790 * np.exp((-1230 - T) * T / (36100 + 360 * T))
    mu_g = 12100 * np.exp((-1233 + T) * T / (9900 + 70 * T))
    a = 0.705 - 0.0017 * T
    b = (4.9 + 0.036 * T) * a ** 2.5
    alpha = 1 - C + a * b * C * (1 - C) / (a * C + b * (1 - C))
    return mu_w ** alpha * mu_g ** (1 - alpha) * 1e-3


def difusion(T, gamma):
    """Einstein: D = kB T / gamma [m²/s]."""
    return KB * T / gamma


def eta_desde_D(D, T, R):
    """eta = kB T / (6 pi R D) [Pa·s]."""
    return KB * T / (6.0 * np.pi * R * D)


def factor_blur_varianza(alpha):
    """S(alpha) de Wong & Halvorsen (Opt. Express 14, 12517, 2006), ec. del encargo.

    <x²>_medida = <x²> S(alpha),  S(alpha) = 2/alpha - 2/alpha² (1 - e^-alpha),  alpha = W / tau_c.
    Es la varianza de la posición PROMEDIADA durante la exposición W de un proceso de
    Ornstein-Uhlenbeck, relativa a la varianza instantánea. S(0) = 1, S decrece con alpha.
    alpha: adimensional (escalar o array). Para alpha < 1e-4 usa la serie 1 - alpha/3 + alpha²/12.
    """
    a = np.asarray(alpha, float)
    chico = np.abs(a) < 1e-4
    ac = np.where(chico, 1.0, a)
    s = 2.0 * (ac + np.expm1(-ac)) / ac ** 2
    return np.where(chico, 1 - a / 3 + a * a / 12, s)


def factor_blur_covarianza(alpha):
    """B(alpha) = 2 (cosh alpha - 1) / alpha² = (sinh(alpha/2) / (alpha/2))².

    Covarianza entre dos posiciones promediadas durante W de un OU separadas tau >= W:
    Cov(x̄(t), x̄(t+tau)) = <x²> B(alpha) e^(-tau/tau_c)   (sale de integrar e^(-|s-s'|/tau_c)
    sobre los dos intervalos de exposición; con alpha -> 0, B -> 1).
    Lo usan la autocorrelación, el espectro y el MSD atrapado de `trampa.py`.
    """
    a = np.asarray(alpha, float) / 2.0
    chico = np.abs(a) < 1e-6
    ac = np.where(chico, 1.0, a)
    return np.where(chico, 1 + a * a / 3, (np.sinh(ac) / ac) ** 2)
