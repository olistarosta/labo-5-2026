"""Tasa de muones vs. ángulo polar (barrido angular) y ajuste R(θ) = R0 cos^n θ.

Para cada carpeta grupo_7_jueves/medidas/polar_<θ>_deg[_vN]:
  N  = eventos seleccionados (coincidencia triple)
  T  = tiempo de medición, desde el inicio de la corrida (nombre run_YYYYMMDD_HHMMSS)
       hasta max(último evento, inicio del último archivo). drs_exam rota de archivo
       cada 300 s, así que T queda determinado a menos de ~1 min (<4 %, despreciable
       frente al error de Poisson, 15-40 %).
  R  = N / T,  σ_R = √N / T  (sólo para graficar)

Ajuste: máxima verosimilitud de Poisson sobre los conteos, μ_i = R0 cos^n(θ_i) T_i.
Es una regresión de Poisson con link log (log μ = log R0 + n log cos θ + log T), así
que se resuelve exacto por Newton y la covarianza es la inversa de la información de
Fisher. No se usa χ² con σ = √N: con pocos conteos subestima los puntos bajos y
sesga n hacia arriba.

Uso:  python fig_barrido_angular.py   -> imprime el ajuste y guarda barrido_angular.{pdf,png,svg}
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2

AQUI = Path(__file__).resolve().parent
DIR_MEDIDAS = AQUI.parent / "grupo_7_jueves" / "medidas"
SEG_POR_ARCHIVO = 300            # MAX_SECS_PER_FILE en drs_exam.cpp
EXCLUIDOS = {15.0}               # ángulos fuera del ajuste (15° norte, las dos mediciones)

_RE_EVENTO = re.compile(r"^Event #(\d+)\s+ts=(\S+ \S+)")
_RE_CARPETA = re.compile(r"polar_(-?[\d.]+)_deg(?:_v(\d+))?$")
_RE_ARCHIVO = re.compile(r"(run_(\d{8}_\d{6}))_idx(\d+)_selected\.txt$")


# ----------------------------------------------------------------- datos
def _leer_corrida(archivos):
    ts = set()
    for a in archivos:
        with open(a) as f:
            for linea in f:
                m = _RE_EVENTO.match(linea)
                if m:
                    ts.add(pd.Timestamp(m.group(2)))
    return sorted(ts)


def cargar(dir_medidas=DIR_MEDIDAS):
    """Una fila por carpeta polar_*: ángulo, N, T [min], tasa y su error [1/min]."""
    filas = []
    for carpeta in sorted(dir_medidas.glob("polar_*_deg*")):
        mc = _RE_CARPETA.match(carpeta.name)
        if not mc:
            continue
        corridas = {}
        for a in carpeta.glob("*_selected.txt"):
            ma = _RE_ARCHIVO.match(a.name)
            corridas.setdefault(ma.group(2), []).append(a)
        N, T, dia = 0, 0.0, None
        for inicio, archivos in sorted(corridas.items()):
            ts = _leer_corrida(archivos)
            if not ts:          # corrida vacía (arranque fallido de polar_0_deg_v2)
                continue
            t0 = pd.to_datetime(inicio, format="%Y%m%d_%H%M%S")
            t_ult = (ts[-1] - t0).total_seconds()
            N += len(ts)
            T += max(t_ult, SEG_POR_ARCHIVO * (len(archivos) - 1))
            dia = t0.strftime("%d/%m")
        theta = float(mc.group(1))
        filas.append(dict(angulo_deg=theta, carpeta=carpeta.name, dia=dia, N=N,
                          T_min=T / 60, tasa=N / (T / 60), err_tasa=np.sqrt(N) / (T / 60),
                          en_ajuste=theta not in EXCLUIDOS))
    return pd.DataFrame(filas)


# ----------------------------------------------------------------- ajustes
def _bondad(N, mu):
    pearson = np.sum((N - mu) ** 2 / mu)
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = 2 * np.sum(mu - N + np.where(N > 0, N * np.log(N / mu), 0))
    return pearson, dev


def ajuste_poisson(df, n_fijo=None):
    """ML de Poisson de R0 cos^n θ. Devuelve dict con R0 [1/min], n, errores, cov de
    (log R0, n), χ² de Pearson, deviance, grados de libertad y p-valor."""
    d = df[df["en_ajuste"]]
    N, T = d["N"].to_numpy(float), d["T_min"].to_numpy()
    x = np.log(np.cos(np.deg2rad(d["angulo_deg"].to_numpy())))
    X = np.c_[np.ones_like(x), x] if n_fijo is None else np.ones((len(x), 1))
    off = np.log(T) + (0 if n_fijo is None else n_fijo * x)
    b = np.zeros(X.shape[1]); b[0] = np.log(N.sum() / T.sum())
    for _ in range(50):                                    # Newton (= IRLS)
        mu = np.exp(X @ b + off)
        paso = np.linalg.solve(X.T @ (mu[:, None] * X), X.T @ (N - mu))
        b += paso
        if np.max(np.abs(paso)) < 1e-12:
            break
    mu = np.exp(X @ b + off)
    cov = np.linalg.inv(X.T @ (mu[:, None] * X))
    if n_fijo is not None:                                 # completar con n sin error
        b = np.r_[b, n_fijo]; cov = np.pad(cov, ((0, 1), (0, 1)))
    pearson, dev = _bondad(N, mu)
    dof = len(N) - X.shape[1]
    R0 = np.exp(b[0])
    return dict(R0=R0, err_R0=R0 * np.sqrt(cov[0, 0]), n=b[1], err_n=np.sqrt(cov[1, 1]),
                cov=cov, beta=b, chi2=pearson, dev=dev, dof=dof, p=chi2.sf(pearson, dof),
                corr=cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]) if n_fijo is None else 0.0)


def ajuste_desplazado(df):
    """Igual pero con el cero desplazado: R0 cos^n(θ - θ0). Sirve para medir la
    asimetría norte/sur. Covarianza por hessiano numérico."""
    d = df[df["en_ajuste"]]
    N, T, th = d["N"].to_numpy(float), d["T_min"].to_numpy(), d["angulo_deg"].to_numpy()

    def mu(p):
        return p[0] * np.cos(np.deg2rad(th - p[2])) ** p[1] * T

    def nll(p):
        m = mu(p)
        return np.sum(m - N * np.log(m))

    p = minimize(nll, [N.sum() / T.sum(), 2, 0], method="Nelder-Mead",
                 options=dict(xatol=1e-10, fatol=1e-12, maxiter=20000)).x
    h = 1e-4 * np.maximum(1, np.abs(p))
    H = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            ei, ej = np.eye(3)[i] * h[i], np.eye(3)[j] * h[j]
            H[i, j] = (nll(p + ei + ej) - nll(p + ei - ej) - nll(p - ei + ej)
                       + nll(p - ei - ej)) / (4 * h[i] * h[j])
    err = np.sqrt(np.diag(np.linalg.inv(H)))
    pearson, dev = _bondad(N, mu(p))
    dof = len(N) - 3
    return dict(R0=p[0], err_R0=err[0], n=p[1], err_n=err[1], theta0=p[2], err_theta0=err[2],
                chi2=pearson, dev=dev, dof=dof, p=chi2.sf(pearson, dof))


def curva(aj, theta_deg):
    """Modelo y banda de 1σ (propagando la covarianza de log R0 y n)."""
    x = np.log(np.cos(np.deg2rad(theta_deg)))
    logR = aj["beta"][0] + aj["beta"][1] * x
    J = np.c_[np.ones_like(x), x]
    s = np.sqrt(np.einsum("ij,jk,ik->i", J, aj["cov"], J))
    return np.exp(logR), np.exp(logR - s), np.exp(logR + s)


# ----------------------------------------------------------------- figura
def figura(df, aj, archivo=None):
    """Mismo estilo que las figuras del informe de láser (seaborn, serif 14 pt)."""
    import matplotlib
    if archivo is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("seaborn-v0_8")
    fs = 14
    plt.rcParams.update({"font.family": "serif", "font.size": fs, "axes.labelsize": fs,
                         "xtick.labelsize": fs, "ytick.labelsize": fs, "legend.fontsize": fs,
                         "mathtext.default": "regular", "mathtext.fontset": "cm",
                         "legend.frameon": True})
    SUR, NORTE, AJUSTE = "#2A9D8F", "#E76F51", "#264653"

    # Serie sur: θ < 0 y el 0° del 29/10 (se midieron juntos); serie norte: θ > 0 y el
    # 0° del 22/10. Los dos 0° se corren 0.8° hacia su lado para que no se tapen.
    sur = (df["angulo_deg"] < 0) | ((df["angulo_deg"] == 0) & (df["dia"] == "29/10"))
    x = df["angulo_deg"] + np.where(df["angulo_deg"] == 0, np.where(sur, -0.8, 0.8), 0)
    rep15 = df.loc[~df["en_ajuste"]].index                   # las dos de 15°: ±0.8°
    x[rep15] += np.linspace(-0.8, 0.8, len(rep15))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    th = np.linspace(-60, 60, 400)
    ax.plot(th, curva(aj, th)[0], color=AJUSTE, lw=2, zorder=2,
            label=r"Ajuste $R_0\,\cos^n\theta$")

    series = [(sur & df["en_ajuste"], SUR, "D", "full", "Hacia el sur"),
              (~sur & df["en_ajuste"], NORTE, "o", "full", "Hacia el norte"),
              (~df["en_ajuste"], NORTE, "o", "hueco", "15° N (fuera del ajuste)")]
    for m, c, mk, relleno, label in series:
        ax.errorbar(x[m], df.loc[m, "tasa"], yerr=df.loc[m, "err_tasa"], fmt=mk,
                    ms=7 if mk == "o" else 6.5, capsize=3, elinewidth=1.6, capthick=1.6,
                    color=c, mfc=c if relleno == "full" else "white", mec=c, mew=1.6,
                    zorder=3, label=label)

    ticks = np.arange(-60, 61, 15)
    ax.set_xticks(ticks, [f"{abs(t)}° {'S' if t < 0 else 'N'}" if t else "0°" for t in ticks])
    ax.set_xlabel(r"Ángulo polar $\theta$")
    ax.set_ylabel("Tasa [eventos/min]")
    ax.set_xlim(-62, 62)
    ax.set_ylim(0, 1.6)
    hl = dict(zip(*reversed(ax.get_legend_handles_labels())))
    orden = ["Hacia el sur", "Hacia el norte", "15° N (fuera del ajuste)",
             r"Ajuste $R_0\,\cos^n\theta$"]
    ax.legend([hl[k] for k in orden], orden, loc="lower center", ncol=2, fontsize=12,
              facecolor="white", framealpha=0.9, columnspacing=1.2)
    fig.tight_layout()
    if archivo is not None:
        for ext in ("pdf", "png", "svg"):
            fig.savefig(f"{archivo}.{ext}", bbox_inches="tight", dpi=300)
    return fig


def resumen(df):
    """Tabla de variantes del ajuste, para decidir qué reportar."""
    filas = []
    todo = df.assign(en_ajuste=True)
    for nombre, d, f in [("sin 15° N", df, ajuste_poisson),
                         ("sin 15° N, n = 2 fijo", df, lambda d: ajuste_poisson(d, n_fijo=2)),
                         ("sin 15° N, cero libre θ0", df, ajuste_desplazado),
                         ("todos", todo, ajuste_poisson),
                         ("todos, cero libre θ0", todo, ajuste_desplazado)]:
        a = f(d)
        filas.append({"ajuste": nombre, "R0 [1/min]": f"{a['R0']:.2f} ± {a['err_R0']:.2f}",
                      "n": f"{a['n']:.2f} ± {a['err_n']:.2f}" if a["err_n"] > 0 else "2 (fijo)",
                      "θ0 [°]": f"{a['theta0']:.1f} ± {a['err_theta0']:.1f}" if "theta0" in a else "—",
                      "χ²/ν": f"{a['chi2']:.1f}/{a['dof']}", "p": round(a["p"], 2)})
    return pd.DataFrame(filas)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(AQUI))
    df = cargar()
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    aj = ajuste_poisson(df)
    print(f"\nR0 = {aj['R0']:.3f} ± {aj['err_R0']:.3f} 1/min   n = {aj['n']:.2f} ± {aj['err_n']:.2f}"
          f"   corr = {aj['corr']:.2f}   χ²/ν = {aj['chi2']:.1f}/{aj['dof']}  (p = {aj['p']:.2f})\n")
    print(resumen(df).to_string(index=False))
    figura(df, aj, AQUI / "barrido_angular")
