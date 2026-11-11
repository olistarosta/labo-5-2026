"""Reporte de decaimientos: cuanto saca cada filtro y como cambia la vida media.

Aplica los filtros de `decaimientos.py` de a uno, en orden, y en cada paso ajusta el
tiempo entre pulsos con una exponencial y con dos. Despues repite moviendo los dos cortes
que mas pesan (amplitud minima del pulso tardio y dt minimo). Imprime las tablas y guarda
las figuras en `muones/salida_decaimientos/reporte_*.pdf` (y .png).

Uso (desde muones/analisis):  python reporte_filtros.py
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import decaimientos as dec

SALIDA = Path(__file__).resolve().parent.parent / "salida_decaimientos"
RANGO = (100.0, 1400.0)            # ventana de dt del ajuste [ns]
ORDEN = ["traza limpia", "entrada en u1", "hay pulso tardio", "amplitud tardio", "dt en rango",
         "ancho tardio", "sin pulsos extra", "sin pulsos previos"]
OPCIONALES = ["dos pulsos en un canal", "tardio donde freno"]

# estilo de ESTILO_FIGURAS.md
plt.style.use("seaborn-v0_8")
fs = 14
plt.rcParams.update({"font.family": "serif", "font.size": fs, "axes.labelsize": fs, "xtick.labelsize": fs,
                     "ytick.labelsize": fs, "legend.fontsize": fs, "mathtext.default": "regular",
                     "mathtext.fontset": "cm", "legend.frameon": True})
C1, C2, CA = "#2A9D8F", "#E76F51", "#264653"
BARRAS = dict(ms=7, capsize=3, elinewidth=1.6, capthick=1.6, mew=1.6, zorder=3)


def guardar(fig, nombre):
    fig.tight_layout()
    fig.savefig(SALIDA / f"{nombre}.pdf", bbox_inches="tight")
    fig.savefig(SALIDA / f"{nombre}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def ajustes(dt, rango, libre=False):
    """Una exponencial, dos con la rapida fija en 206 ns y (si `libre`) dos con todo libre."""
    r = {"una": dec.ajustar_vida_media(dt, rango, 0.0),
         "fija": dec.ajustar_dos_exponenciales(dt, rango, tau_lenta=None, tau_rapida=dec.TAU_MU_FE)}
    if libre:
        r["libre"] = dec.ajustar_dos_exponenciales(dt, rango, tau_lenta=None, tau_rapida=None)
    return r


def fila(nombre, n, r):
    u, f = r["una"], r["fija"]
    txt = (f"{nombre:26s} {n:7d} | {u['tau']:5.0f} -{u['err_menos']:3.0f} +{u['err_mas']:3.0f} | "
           f"{f['tau_lenta']:5.0f}  {100 * f['f']:4.1f} [{100 * f['f_min']:4.1f},{100 * f['f_max']:5.1f}]  {f['sigmas']:3.1f}")
    if "libre" in r:
        l = r["libre"]
        txt += (f" | {l['tau_lenta']:6.0f} {l['tau_rapida']:5.0f}  {100 * l['f']:4.1f} "
                f"[{100 * l['f_min']:4.1f},{100 * l['f_max']:5.1f}]  {l['sigmas'] ** 2:4.1f}")
    return txt


CAB = (f"{'':26s} {'eventos':>7s} | {'una exp: tau [ns]':17s} | dos exp, rapida = 206 ns fija: "
       "tau lenta, % rapida [68 %], mejora en sigmas")
CAB_LIBRE = " | dos exp libres: tau lenta, tau rapida, % rapida [68 %], 2*Delta(lnL)"


def main():
    ev = dec.cargar_indice()
    tab = dec.tabla_eventos(ev)
    vetos = dec.mascara_vetos(tab, dec.cargar_vetos(SALIDA / "vetos_manuales.csv"))
    dt = tab["dt"].to_numpy()
    filtros = dec.armar_filtros(tab, dt_rango=RANGO)

    # ---------------------------------------------------------------- 1. filtros de a uno
    print(f"\n=== 1. Filtros en orden (ventana de ajuste {RANGO[0]:.0f}-{RANGO[1]:.0f} ns) ===")
    print(f"{'filtro':26s} {'pasan solo':>10s} {'quedan':>8s} {'saca':>8s} {'% que saca':>10s}")
    m = vetos.copy()
    pasos, quedan = ["todos"], [int(m.sum())]
    print(f"{'todos los eventos':26s} {'':>10s} {m.sum():8d}")
    for nombre in ORDEN + OPCIONALES:
        antes = m.sum()
        m = m & filtros[nombre]
        pasos.append(nombre)
        quedan.append(int(m.sum()))
        print(f"{nombre:26s} {filtros[nombre].sum():10d} {m.sum():8d} {antes - m.sum():8d} {100 * (antes - m.sum()) / antes:9.1f}%")

    print("\n" + CAB + CAB_LIBRE)
    m = vetos.copy()
    res_pasos = {}
    for nombre in ORDEN + OPCIONALES:
        m = m & filtros[nombre]
        if ORDEN.index("hay pulso tardio") <= (ORDEN + OPCIONALES).index(nombre):
            en_rango = m & (dt >= RANGO[0]) & (dt <= RANGO[1])
            res_pasos[nombre] = (int(en_rango.sum()), ajustes(dt[m], RANGO, libre=True))
            print(fila("+ " + nombre, *res_pasos[nombre]), flush=True)

    # ---------------------------------------------------------------- 2. mover los dos cortes principales
    def seleccion(**cambios):
        f = dec.armar_filtros(tab, **{"dt_rango": RANGO, **cambios})
        return dec.combinar(f, ORDEN) & vetos

    print("\n=== 2. Amplitud minima del pulso tardio (resto de los filtros igual) ===\n" + CAB)
    res_amp = {}
    for a in (30, 40, 50, 60, 70, 80, 100):
        s = seleccion(a_tardio_min=a)
        res_amp[a] = (int(s.sum()), ajustes(dt[s], RANGO))
        print(fila(f"amplitud >= {a} mV", *res_amp[a]), flush=True)

    print("\n=== 3. dt minimo (resto de los filtros igual; la ventana termina en 1400 ns) ===\n" + CAB)
    res_dt = {}
    for a in (60, 80, 100, 150, 200, 300):
        s = seleccion(dt_rango=(a, RANGO[1]))
        res_dt[a] = (int(s.sum()), ajustes(dt[s], (a, RANGO[1])))
        print(fila(f"dt >= {a} ns", *res_dt[a]), flush=True)

    print("\n=== 4. Seleccion final, tres ventanas ===\n" + CAB + CAB_LIBRE)
    res_final = {}
    for a in (60, 100, 200):
        s = seleccion(dt_rango=(a, RANGO[1]))
        res_final[a] = (int(s.sum()), ajustes(dt[s], (a, RANGO[1]), libre=True), dt[s])
        print(fila(f"dt en {a}-1400 ns", *res_final[a][:2]), flush=True)

    # ---------------------------------------------------------------- 5. que pasa por debajo de 200 ns
    # a) el exceso de candidatos en 60-200 ns, segun el pulso de entrada del canal donde aparece el tardio.
    #    Un pospulso o la cola del pulso de entrada solo pueden estar en un canal que tuvo entrada, y
    #    crecen con su amplitud.
    s = seleccion(dt_rango=(60, RANGO[1]))
    canal = np.maximum(tab["canal_tardio"].to_numpy() - 1, 0)
    a_ent = np.stack([tab[f"a_ent_{c}"].to_numpy() for c in dec.CANALES], axis=1)[np.arange(len(tab)), canal]
    tau = dec.TAU_MU
    escala = (np.exp(-60 / tau) - np.exp(-200 / tau)) / (np.exp(-300 / tau) - np.exp(-1400 / tau))
    print("\n=== 5. Exceso de candidatos en 60-200 ns (esperado: exponencial de 2197 ns normalizada en 300-1400 ns) ===")
    for nombre, m in (("tardio en un canal sin entrada", a_ent == 0), ("entrada < 150 mV", (a_ent > 0) & (a_ent < 150)),
                      ("entrada 150-300 mV", (a_ent >= 150) & (a_ent < 300)), ("entrada > 300 mV", a_ent >= 300),
                      ("todos", np.ones(len(tab), bool))):
        x = dt[s & m]
        obs, ref = ((x >= 60) & (x < 200)).sum(), ((x >= 300) & (x <= 1400)).sum()
        print(f"{nombre:32s} observado {obs:4d}   esperado {ref * escala:6.1f}   exceso {obs - ref * escala:+6.1f} ± {np.sqrt(obs):.1f}")

    # a2) cuantos pospulsos hay, usando una muestra de control: si el tardio aparece en un centellador
    #     que no tuvo pulso de entrada no puede ser un pospulso. La proporcion cortos/largos de esa
    #     muestra dice cuantos cortos reales esperar en la otra; lo que sobre son pospulsos.
    pe, pt = tab["patron_ent"].to_numpy(), tab["patron_tar"].to_numpy()
    CORTO, LARGO = (60, 200), (300, 1400)

    def cuenta(m, r):
        return int((s & m & (dt >= r[0]) & (dt < r[1])).sum())

    print(f"\n=== 5b. Pospulsos por comparacion con un control (cortos = {CORTO[0]}-{CORTO[1]} ns, largos = {LARGO[0]}-{LARGO[1]} ns) ===")
    for nombre, control, prueba in (
            ("frena en u2 o hierro: electron solo en u3 (control) contra electron solo en u2", (pe == "110") & (pt == "001"), (pe == "110") & (pt == "010")),
            ("todos: tardio en canal sin entrada (control) contra canal con entrada", a_ent == 0, a_ent > 0)):
        cc, cl, pc, pl = cuenta(control, CORTO), cuenta(control, LARGO), cuenta(prueba, CORTO), cuenta(prueba, LARGO)
        esperado = cc / cl * pl
        exceso, err = pc - esperado, np.sqrt(pc + esperado ** 2 * (1 / cc + 1 / pl))
        print(f"{nombre}\n   control {cc}/{cl} = {cc / cl:.3f} ± {cc / cl * np.sqrt(1 / cc + 1 / cl):.3f}   prueba {pc}/{pl} = {pc / pl:.3f} ± "
              f"{pc / pl * np.sqrt(1 / pc + 1 / pl):.3f}   (exponencial pura: {escala:.3f})")
        print(f"   cortos esperados en la prueba sin pospulsos: {esperado:.0f}; observados {pc}  ->  pospulsos = {exceso:+.0f} ± {err:.0f}"
              f"   (cota 95 %: menos de {max(exceso, 0) + 1.64 * err:.0f}, el {100 * (max(exceso, 0) + 1.64 * err) / pc:.0f} % de los cortos)")
    bordes_c = np.array([60, 100, 150, 200, 300, 500, 800, 1100, 1400])
    control_prueba = {}
    for clave, m in (("con", a_ent > 0), ("sin", a_ent == 0)):
        n_bin = np.histogram(dt[s & m], bordes_c)[0]
        control_prueba[clave] = (n_bin, cuenta(m, LARGO))

    # a3) viene del hierro? Los mu- que frenan en la placa de hierro tienen vida media ~206 ns (casi todos
    #     son capturados; decae ~9 %). Si es eso, el exceso esta en los que pudieron frenar en el hierro
    #     (entrada 110) y no en los que frenaron en u3 (entrada 111, solo plastico), y su tamano se puede
    #     predecir sin parametros libres.
    A_FE, B_FE = 60.0, RANGO[1]
    ventana = lambda t_, a_=A_FE, b_=B_FE: np.exp(-a_ / t_) - np.exp(-b_ / t_)
    F_MAS, F_DECAE_FE = 0.56, dec.TAU_MU_FE / tau * 0.975
    rapidos_por_lento = (1 - F_MAS) / F_MAS * F_DECAE_FE * ventana(dec.TAU_MU_FE) / ventana(tau)
    f_hierro = rapidos_por_lento / (1 + rapidos_por_lento)          # fraccion rapida en eventos del hierro
    clases_fe = {"frena en u3 (solo plastico)": pe == "111",
                 "frena en u2 o en el hierro": pe == "110",
                 "   electron solo en u3 (hierro seguro)": (pe == "110") & (pt == "001"),
                 "   electron solo en u2": (pe == "110") & (pt == "010")}
    n_abajo = int((s & clases_fe["   electron solo en u3 (hierro seguro)"]).sum())
    n_110 = int((s & (pe == "110")).sum())
    f_110 = f_hierro * 2 * n_abajo / n_110            # del hierro salen tantos electrones hacia arriba como hacia abajo
    prediccion = {"frena en u3 (solo plastico)": 0.0, "frena en u2 o en el hierro": f_110,
                  "   electron solo en u3 (hierro seguro)": f_hierro, "   electron solo en u2": np.nan}
    print(f"\n=== 5c. Captura de mu- en el hierro (ventana {A_FE:.0f}-{B_FE:.0f} ns) ===")
    print(f"{'clase':40s} {'n':>5s} {'exceso en 60-200 ns':>20s} {'rapida 206 ns [%]':>24s} {'prediccion [%]':>15s}")
    res_fe = {}
    for nombre, m in clases_fe.items():
        x = dt[s & m]
        obs, ref = ((x >= 60) & (x < 200)).sum(), ((x >= 300) & (x <= 1400)).sum()
        r = dec.ajustar_dos_exponenciales(x, (A_FE, B_FE), tau_lenta=tau, tau_rapida=dec.TAU_MU_FE)
        res_fe[nombre] = r
        print(f"{nombre:40s} {len(x):5d} {obs - ref * escala:+10.1f} ± {np.sqrt(obs):4.1f} "
              f"{100 * r['f']:8.1f} [{100 * r['f_min']:.1f}, {100 * r['f_max']:.1f}] ({r['sigmas']:.1f} sigma) {100 * prediccion[nombre]:12.1f}")
    print(f"prediccion de exceso en 60-200 ns para 'frena en u2 o en el hierro': "
          f"{f_110 * n_110 * ventana(dec.TAU_MU_FE, 60, 200) / ventana(dec.TAU_MU_FE):.0f} eventos")
    clases_dt = {k: (np.histogram(dt[s & m], bordes_c)[0], cuenta(m, LARGO)) for k, m in (("110", pe == "110"), ("111", pe == "111"))}

    # b) directo sobre las trazas, sin el buscador de pulsos: promedio de muchas trazas alineadas en el
    #    pulso de entrada, y cuantas bajan de -20 mV entre 60 y 200 ns despues
    limpio = dec.combinar(filtros, ["traza limpia", "entrada en u1"]) & (tab["patron_ent"].to_numpy() == "111")
    elegidos = np.sort(np.random.default_rng(0).choice(np.flatnonzero(limpio), 12000, replace=False))
    paso, antes, despues = 1.66, 60, 300
    suma, n_trazas = np.zeros((antes + despues, 3)), 0
    minimos, entradas = {"60-200": [], "600-740": []}, []
    t0 = tab["t0"].to_numpy()
    for i in elegidos:
        t, v = ev.onda(i)
        j0 = int(round(t0[i] / paso))
        if len(t) != 1024 or j0 < antes:
            continue
        x, xs, _ = dec.senales(v)
        suma += x[j0 - antes:j0 + despues]
        n_trazas += 1
        entradas.append([tab[f"a_ent_{c}"].iloc[i] for c in dec.CANALES])
        for clave in minimos:
            lo, hi = (int(float(q) / paso) for q in clave.split("-"))
            minimos[clave].append(xs[j0 + lo:j0 + hi].min(axis=0))
    media = suma / n_trazas
    t_media = (np.arange(antes + despues) - antes) * paso
    entradas = np.array(entradas)
    minimos = {k: np.array(v) for k, v in minimos.items()}
    print(f"\n=== 6. Cola del pulso de entrada ({n_trazas} trazas de muones que atraviesan) ===")
    print("traza media despues del pulso de entrada [mV], u1 u2 u3:")
    for lo, hi in ((30, 50), (50, 80), (80, 120), (120, 200), (200, 300), (300, 450)):
        print(f"   {lo:3d}-{hi:3d} ns: ", media[(t_media >= lo) & (t_media < hi)].mean(axis=0).round(2))
    print("trazas que bajan de -20 mV (suavizado) entre 60 y 200 ns despues de la entrada, segun la amplitud de entrada del canal:")
    for lo, hi in ((60, 150), (150, 250), (250, 600)):
        txt = []
        for c in range(3):
            m = (entradas[:, c] >= lo) & (entradas[:, c] < hi)
            txt.append(f"u{c + 1}: {100 * (minimos['60-200'][m, c] < -20).mean():5.2f}% (control 600-740 ns: "
                       f"{100 * (minimos['600-740'][m, c] < -20).mean():4.2f}%)")
        print(f"   entrada {lo:3d}-{hi:3d} mV:  " + "   ".join(txt))

    # ---------------------------------------------------------------- figuras
    # a) cuantos eventos quedan despues de cada filtro
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(pasos))
    ax.plot(x[:len(ORDEN) + 1], quedan[:len(ORDEN) + 1], "o", color=C1, ms=8, zorder=3, label="Filtros de la selección")
    ax.plot(x[len(ORDEN) + 1:], quedan[len(ORDEN) + 1:], "o", color=C1, mfc="white", mew=1.6, ms=8, zorder=3,
            label="Filtros opcionales de topología")
    ax.set_yscale("log")
    ax.set_xticks(x, [p.capitalize() for p in pasos], rotation=40, ha="right", fontsize=11)
    ax.set_ylabel("Eventos que quedan")
    ax.legend(loc="upper right", facecolor="white", framealpha=0.9, fontsize=12)
    guardar(fig, "reporte_filtros_eventos")

    # b) tau despues de cada filtro
    def tau_y_error(resultados):
        t = np.array([r["una"]["tau"] for _, r in resultados])
        e = np.array([[r["una"]["err_menos"], r["una"]["err_mas"]] for _, r in resultados]).T
        return t, e

    nombres = list(res_pasos)
    t, e = tau_y_error(res_pasos.values())
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.axhline(dec.TAU_MU, color=CA, lw=2, zorder=2, label="Muón libre")
    k = len(nombres) - len(OPCIONALES)
    ax.errorbar(np.arange(k), t[:k], yerr=e[:, :k], fmt="o", color=C1, mfc=C1, mec=C1, label="Filtros de la selección", **BARRAS)
    ax.errorbar(np.arange(k, len(nombres)), t[k:], yerr=e[:, k:], fmt="o", color=C1, mfc="white", mec=C1,
                label="Filtros opcionales de topología", **BARRAS)
    ax.set_xticks(np.arange(len(nombres)), [f"+ {n}\n({res_pasos[n][0]})" for n in nombres], rotation=40, ha="right", fontsize=10)
    ax.set_ylabel(r"Vida media $\tau$ [ns]")
    ax.set_ylim(0, 4000)
    ax.legend(loc="lower right", facecolor="white", framealpha=0.9, fontsize=12)
    guardar(fig, "reporte_filtros_tau")

    # c) y d) tau al mover los dos cortes principales
    for res, etiqueta, nombre in ((res_amp, "Amplitud mínima del pulso tardío [mV]", "reporte_tau_vs_amplitud"),
                                  (res_dt, r"$\Delta t$ mínimo [ns]", "reporte_tau_vs_dtmin")):
        t, e = tau_y_error(res.values())
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.axhline(dec.TAU_MU, color=CA, lw=2, zorder=2, label="Muón libre")
        ax.errorbar(list(res), t, yerr=e, fmt="o", color=C1, mfc=C1, mec=C1, label="Ajuste de una exponencial", **BARRAS)
        ax.set_xlabel(etiqueta)
        ax.set_ylabel(r"Vida media $\tau$ [ns]")
        ax.set_ylim(0, 5000)
        ax.legend(loc="upper left", facecolor="white", framealpha=0.9, fontsize=12)
        guardar(fig, nombre)

    # e) histograma de dt con los dos ajustes (ventana desde 60 ns, donde pesa la componente rapida)
    n, r, x = res_final[60]
    a, b = 60.0, RANGO[1]
    x = x[(x >= a) & (x <= b)]
    bordes = np.linspace(a, b, 21)
    cuentas = np.histogram(x, bordes)[0]
    ancho = bordes[1] - bordes[0]
    tt = np.linspace(a, b, 500)
    una = len(x) * ancho * dec._exp_norm(tt, r["una"]["tau"], a, b)
    f = r["fija"]
    dos = len(x) * ancho * (f["f"] * dec._exp_norm(tt, f["tau_rapida"], a, b) + (1 - f["f"]) * dec._exp_norm(tt, f["tau_lenta"], a, b))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(tt / 1000, una, color=CA, lw=2, zorder=2, label="Una exponencial")
    ax.plot(tt / 1000, dos, color=C2, lw=2, ls="--", zorder=2, label="Dos exponenciales")
    ax.errorbar(0.5 * (bordes[1:] + bordes[:-1]) / 1000, cuentas, yerr=np.sqrt(cuentas), fmt="o", color=C1, mfc=C1, mec=C1,
                label="Decaimientos", **BARRAS)
    ax.set_xlabel(r"Tiempo entre pulsos $\Delta t$ [$\mu$s]")
    ax.set_ylabel(f"Eventos / {ancho:.0f} ns")
    ax.set_ylim(0, None)
    ax.legend(loc="upper right", facecolor="white", framealpha=0.9)
    guardar(fig, "reporte_ajuste_una_y_dos_exp")

    # f) traza media despues del pulso de entrada: la cola que llega hasta ~200 ns
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for c, (color, estilo) in enumerate(zip(("#4C72B0", C2, C1), ("-", "--", "-."))):
        ax.plot(t_media, media[:, c], color=color, ls=estilo, lw=2, label=f"Centellador {c + 1}")
    ax.set_xlim(0, 450)
    ax.set_ylim(-8, 2)
    ax.set_xlabel("Tiempo desde el pulso de entrada [ns]")
    ax.set_ylabel("Tensión media [mV]")
    ax.legend(loc="lower right", facecolor="white", framealpha=0.9)
    guardar(fig, "reporte_cola_pulso_entrada")

    # g) dt del tardio segun si su centellador tuvo o no pulso de entrada (control de pospulsos)
    centro = 0.5 * (bordes_c[1:] + bordes_c[:-1])
    por_100 = np.diff(bordes_c) / 100
    tt = np.linspace(bordes_c[0], bordes_c[-1], 300)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(tt, 100 * 100 * np.exp(-tt / tau) / (tau * (np.exp(-LARGO[0] / tau) - np.exp(-LARGO[1] / tau))), color=CA, lw=2,
            zorder=2, label="Exponencial de 2197 ns")
    for clave, color, marca, corrimiento, etiqueta in (("con", C1, "o", -8, "Electrón en un centellador con pulso de entrada"),
                                                       ("sin", C2, "D", +8, "Electrón en un centellador sin pulso de entrada")):
        n_bin, n_largo = control_prueba[clave]
        ax.errorbar(centro + corrimiento, 100 * n_bin / por_100 / n_largo, yerr=100 * np.sqrt(n_bin) / por_100 / n_largo, fmt=marca,
                    color=color, mfc=color, mec=color, label=etiqueta, **BARRAS)
    ax.set_xlabel(r"Tiempo entre pulsos $\Delta t$ [ns]")
    ax.set_ylabel("Eventos cada 100 ns [% de los de 300 a 1400 ns]")
    ax.yaxis.label.set_size(12)
    ax.set_ylim(0, None)
    ax.legend(loc="lower left", facecolor="white", framealpha=0.9, fontsize=11)
    guardar(fig, "reporte_pospulsos_control")

    # h) dt segun donde freno el muon, con la prediccion de la captura en hierro
    def curva(f):
        "Eventos cada 100 ns, en % de los de 300-1400 ns, para una fraccion f de componente de 206 ns."
        forma = lambda x: f * dec._exp_norm(x, dec.TAU_MU_FE, A_FE, B_FE) + (1 - f) * dec._exp_norm(x, tau, A_FE, B_FE)
        en_largos = f * ventana(dec.TAU_MU_FE, *LARGO) / ventana(dec.TAU_MU_FE) + (1 - f) * ventana(tau, *LARGO) / ventana(tau)
        return 100 * 100 * forma(tt) / en_largos

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(tt, curva(0.0), color=CA, lw=2, zorder=2, label="Exponencial de 2197 ns")
    ax.plot(tt, curva(f_110), color=CA, lw=2, ls="--", zorder=2, label=r"Con los $\mu^-$ capturados en el hierro")
    for clave, color, marca, corrimiento, etiqueta in (("110", C1, "o", -8, "Frenan en el centellador 2 o en el hierro"),
                                                       ("111", C2, "D", +8, "Frenan en el centellador 3 (solo plástico)")):
        n_bin, n_largo = clases_dt[clave]
        ax.errorbar(centro + corrimiento, 100 * n_bin / por_100 / n_largo, yerr=100 * np.sqrt(n_bin) / por_100 / n_largo, fmt=marca,
                    color=color, mfc=color, mec=color, label=etiqueta, **BARRAS)
    ax.set_xlabel(r"Tiempo entre pulsos $\Delta t$ [ns]")
    ax.set_ylabel("Eventos cada 100 ns [% de los de 300 a 1400 ns]")
    ax.yaxis.label.set_size(12)
    ax.set_ylim(0, None)
    ax.legend(loc="lower left", facecolor="white", framealpha=0.9, fontsize=11)
    guardar(fig, "reporte_exceso_hierro")
    print("\nfiguras guardadas en", SALIDA)


if __name__ == "__main__":
    main()
