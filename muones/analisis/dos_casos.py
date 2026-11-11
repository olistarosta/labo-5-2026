"""Ajuste de una y de dos exponenciales en dos grupos de decaimientos.

  Grupo 1, frena en el hierro: pulso de entrada en los centelladores 1 y 2 y no en el 3.
  Grupo 2, cruza el hierro:    pulso de entrada en los tres; el muon frena en el plastico
                               del centellador 3.

Para cada grupo se ajusta el tiempo entre pulsos, por maxima verosimilitud sin binear, con
  una exponencial:   exp(-dt/tau)
  dos exponenciales: f*exp(-dt/tau_rapida) + (1-f)*exp(-dt/tau_lenta), todo libre
y se informan los parametros con su intervalo de 68 % y la bondad de cada ajuste
(chi2 sobre el histograma, mejora de la verosimilitud y su p-valor por simulacion).

Como script (desde muones/analisis):  python dos_casos.py
   imprime las tablas y guarda casos_frena_en_hierro y casos_cruza_el_hierro (.pdf y .png)
   en muones/salida_decaimientos/.
Como modulo:  import dos_casos as dc;  grupos = dc.cargar_grupos();  res = dc.analizar(x)
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import chi2

import decaimientos as dec

SALIDA = Path(__file__).resolve().parent.parent / "salida_decaimientos"
FILTROS = ["traza limpia", "entrada en u1", "hay pulso tardio", "amplitud tardio", "dt en rango",
           "ancho tardio", "sin pulsos extra", "sin pulsos previos"]
A, B = 60.0, 1400.0            # ventana de dt [ns]
N_BINS = 20                    # bins del histograma para el chi2
# limites del ajuste de dos exponenciales: tau rapida, tau lenta [ns] y fraccion rapida
LIMITES = [(20.0, 800.0), (800.0, 20000.0), (0.0, 0.95)]
NOMBRES = ("tau_rapida", "tau_lenta", "fraccion")
GRUPOS = {"Grupo 1: frena en el hierro": ("110", "casos_frena_en_hierro"),
          "Grupo 2: cruza el hierro y frena en el centellador 3": ("111", "casos_cruza_el_hierro")}
# puntos de partida del minimizador: la verosimilitud es casi plana y desde un solo punto
# puede frenar antes del minimo, asi que se prueban varios y se toma el de menor -lnL
INICIOS = [(tr, tl, f) for tr in (40.0, 100.0, 200.0, 400.0) for tl in (1500.0, 2200.0, 4000.0) for f in (0.02, 0.1, 0.3)]
INICIOS_POCOS = [(tr, tl, f) for tr in (60.0, 250.0) for tl in (1800.0, 3000.0) for f in (0.03, 0.2)]


# ------------------------------------------------------------------ modelo
def exp_norm(x, tau):
    """Exponencial normalizada en la ventana [A, B]."""
    return np.exp(-x / tau) / (tau * (np.exp(-A / tau) - np.exp(-B / tau)))


def dens(x, p):
    """Densidad del modelo: p = (tau,) o p = (tau_rapida, tau_lenta, fraccion)."""
    if len(p) == 1:
        return exp_norm(x, p[0])
    tr, tl, f = p
    return f * exp_norm(x, tr) + (1 - f) * exp_norm(x, tl)


def menos_lnL(p, x):
    return -np.log(dens(x, np.asarray(p))).sum()


# ------------------------------------------------------------------ ajustes
def ajustar_dos(x, fijos=None, inicios=INICIOS):
    """Minimo de -lnL para dos exponenciales. `fijos` = {indice: valor} deja parametros quietos."""
    mejor = None
    for p0 in inicios:
        p0, lim = list(p0), list(LIMITES)
        for k, v in (fijos or {}).items():
            p0[k], lim[k] = v, (v, v)
        r = minimize(menos_lnL, p0, args=(x,), method="L-BFGS-B", bounds=lim, options={"ftol": 1e-12, "gtol": 1e-9})
        if mejor is None or r.fun < mejor.fun:
            mejor = r
    return mejor


def intervalo(x, k, mejor):
    """Intervalo de 68 % del parametro k: donde el perfil de -lnL sube menos de 1/2."""
    lo, hi = LIMITES[k]
    if k == 2:                                   # fraccion: rejilla fina donde suele estar
        rejilla = np.concatenate([np.linspace(lo, 0.1, 41), np.linspace(0.1, hi, 35)[1:]])
    else:
        rejilla = np.geomspace(lo, hi, 41)
    inicios = INICIOS_POCOS + [tuple(mejor.x)]
    dentro = [v for v in rejilla if ajustar_dos(x, {k: v}, inicios).fun - mejor.fun <= 0.5] + [mejor.x[k]]
    return min(dentro), max(dentro)


def chi2_binado(x, p):
    """chi2 de Pearson del histograma contra el modelo, grados de libertad y p-valor."""
    bordes = np.linspace(A, B, N_BINS + 1)
    obs = np.histogram(x, bordes)[0]
    fina = np.linspace(A, B, 20 * N_BINS + 1)                      # integral del modelo en cada bin
    d = dens(fina, np.asarray(p))
    esp = len(x) * ((d[1:] + d[:-1]) / 2 * np.diff(fina)).reshape(N_BINS, 20).sum(axis=1)
    c2 = float(((obs - esp) ** 2 / esp).sum())
    gl = N_BINS - 1 - len(p)
    return c2, gl, float(chi2.sf(c2, gl))


def simular(n, tau, rng):
    """n tiempos de una exponencial de vida media tau, dentro de la ventana [A, B]."""
    u = rng.uniform(size=n)
    return -tau * np.log(np.exp(-A / tau) - u * (np.exp(-A / tau) - np.exp(-B / tau)))


def p_simulado(n, tau, mejora, n_sim=200):
    """Fraccion de muestras simuladas con UNA sola exponencial en las que agregar la segunda
    mejora la verosimilitud tanto o mas que en los datos. La formula habitual (chi2 con 2
    grados de libertad) no vale aca, porque la fraccion no puede ser negativa."""
    rng = np.random.default_rng(0)
    iguala = 0
    for _ in range(n_sim):
        x = simular(n, tau, rng)
        una = dec.ajustar_vida_media(x, (A, B))
        iguala += 2 * (menos_lnL((una["tau"],), x) - ajustar_dos(x).fun) >= mejora
    return iguala / n_sim


# ------------------------------------------------------------------ analisis de un grupo
def cargar_grupos():
    """{nombre: tiempos entre pulsos [ns] de los decaimientos de ese grupo, dentro de la ventana}."""
    ev = dec.cargar_indice()
    tab = dec.tabla_eventos(ev)
    sel = dec.combinar(dec.armar_filtros(tab, dt_rango=(A, B)), FILTROS)
    dt, entrada = tab["dt"].to_numpy(), tab["patron_ent"].to_numpy()
    grupos = {}
    for nombre, (patron, _) in GRUPOS.items():
        x = dt[sel & (entrada == patron)]
        grupos[nombre] = x[(x >= A) & (x <= B)]
    return grupos


def analizar(x, n_sim=200):
    """Todos los numeros de un grupo. Devuelve un dict."""
    una = dec.ajustar_vida_media(x, (A, B))
    p1 = (una["tau"],)
    dos = ajustar_dos(x)
    esperado = ajustar_dos(x, {0: dec.TAU_MU_FE, 1: dec.TAU_MU})
    mejora = 2 * (menos_lnL(p1, x) - dos.fun)
    return dict(
        n=len(x), p1=p1, p2=dos.x,
        tau=una["tau"], tau_int=(una["tau"] - una["err_menos"], una["tau"] + una["err_mas"]),
        chi2_una=chi2_binado(x, p1), chi2_dos=chi2_binado(x, dos.x),
        intervalos={NOMBRES[k]: intervalo(x, k, dos) for k in range(3)},
        esperado_fraccion=esperado.x[2], esperado_dlnL=esperado.fun - dos.fun,
        mejora=mejora, p_mejora=p_simulado(len(x), una["tau"], mejora, n_sim),
    )


def texto_intervalo(valor, lo, hi, k):
    """El valor y su intervalo, avisando cuando llega al limite permitido (el dato no lo acota)."""
    num = (lambda v: f"{100 * v:.1f} %") if k == 2 else (lambda v: f"{v:.0f} ns")
    izq = "sin cota" if np.isclose(lo, LIMITES[k][0]) else num(lo)
    der = "sin cota" if np.isclose(hi, LIMITES[k][1]) else num(hi)
    return f"{num(valor)} (68 %: {izq} a {der})"


def imprimir(nombre, r):
    c1, g1, pv1 = r["chi2_una"]
    c2, g2, pv2 = r["chi2_dos"]
    print(f"\n{nombre}: {r['n']} eventos")
    print("  UNA exponencial")
    print(f"     tau         = {r['tau']:.0f} ns (68 %: {r['tau_int'][0]:.0f} a {r['tau_int'][1]:.0f})")
    print(f"     chi2/gl     = {c1:.1f}/{g1} = {c1 / g1:.2f}   p-valor = {pv1:.2f}")
    print("  DOS exponenciales")
    for k, etiqueta in enumerate(("tau rapida ", "tau lenta  ", "fraccion   ")):
        print(f"     {etiqueta} = {texto_intervalo(r['p2'][k], *r['intervalos'][NOMBRES[k]], k)}")
    print(f"     chi2/gl     = {c2:.1f}/{g2} = {c2 / g2:.2f}   p-valor = {pv2:.2f}")
    print(f"  DOS exponenciales con las vidas medias esperadas fijas ({dec.TAU_MU_FE:.0f} y {dec.TAU_MU:.0f} ns)")
    print(f"     fraccion    = {100 * r['esperado_fraccion']:.1f} %")
    print(f"     -lnL        = {r['esperado_dlnL']:+.2f} respecto del mejor ajuste de dos exponenciales"
          f"  (hasta +1,15 es compatible al 68 %, hasta +3,00 al 95 %)")
    print("  MEJORA al pasar de una a dos")
    print(f"     chi2 baja {c1 - c2:.1f} con 2 parametros mas (el chi2/gl puede subir igual, porque se divide por menos)")
    print(f"     2*Delta(lnL) = {r['mejora']:.2f}")
    print(f"     p-valor simulado = {r['p_mejora']:.3f}  (fraccion de muestras de UNA sola exponencial "
          f"que mejoran tanto o mas; menor a 0,05 seria significativo)")


def figura(x, r):
    """Histograma de dt con el ajuste de una y de dos exponenciales (estilo de ESTILO_FIGURAS.md)."""
    bordes = np.linspace(A, B, N_BINS + 1)
    cuentas, ancho = np.histogram(x, bordes)[0], bordes[1] - bordes[0]
    t = np.linspace(A, B, 500)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(t / 1000, len(x) * ancho * dens(t, np.asarray(r["p1"])), color="#264653", lw=2, zorder=2, label="Una exponencial")
    ax.plot(t / 1000, len(x) * ancho * dens(t, r["p2"]), color="#E76F51", lw=2, ls="--", zorder=2, label="Dos exponenciales")
    ax.errorbar(0.5 * (bordes[1:] + bordes[:-1]) / 1000, cuentas, yerr=np.sqrt(cuentas), fmt="o", ms=7, capsize=3,
                elinewidth=1.6, capthick=1.6, mew=1.6, color="#2A9D8F", mfc="#2A9D8F", mec="#2A9D8F", zorder=3,
                label="Decaimientos")
    ax.set_xlabel(r"Tiempo entre pulsos $\Delta t$ [$\mu$s]")
    ax.set_ylabel(f"Eventos / {ancho:.0f} ns")
    ax.set_ylim(0, None)
    ax.legend(loc="upper right", facecolor="white", framealpha=0.9)
    fig.tight_layout()
    return fig


def estilo():
    """Configuracion de ESTILO_FIGURAS.md."""
    plt.style.use("seaborn-v0_8")
    fs = 14
    plt.rcParams.update({"font.family": "serif", "font.size": fs, "axes.labelsize": fs, "xtick.labelsize": fs,
                         "ytick.labelsize": fs, "legend.fontsize": fs, "mathtext.default": "regular",
                         "mathtext.fontset": "cm", "legend.frameon": True})


def main():
    plt.switch_backend("Agg")
    estilo()
    print(f"Ventana de dt: {A:.0f} a {B:.0f} ns. chi2 sobre {N_BINS} bins de {(B - A) / N_BINS:.0f} ns.")
    for nombre, x in cargar_grupos().items():
        r = analizar(x)
        imprimir(nombre, r)
        fig = figura(x, r)
        fig.savefig(SALIDA / f"{GRUPOS[nombre][1]}.pdf", bbox_inches="tight")
        fig.savefig(SALIDA / f"{GRUPOS[nombre][1]}.png", bbox_inches="tight", dpi=300)
        plt.close(fig)
    print("\nfiguras guardadas en", SALIDA)


if __name__ == "__main__":
    main()
