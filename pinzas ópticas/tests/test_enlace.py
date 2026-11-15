"""Tests de pinzas.enlace: detectar en cada cuadro + enlazar + postprocesar, sobre sintéticos.

Videos hechos en memoria con sintetico.dibujar + sintetico.camara (cuentas crudas, fondo
despareja), con la verdad de cada bead por cuadro. Pureza de una trayectoria = fracción de sus
filas medidas que caen (< 1 px) sobre la bead verdadera mayoritaria; tiene que ser 1.
"""
import numpy as np
import pytest

from pinzas import enlace as en
from pinzas import sintetico as sint
from pinzas import tracking as tr

R = 6.0                  # px: bead de ~1.1 µm de diámetro a 90 nm/px
NIVEL, CONTRASTE = 600.0, 0.35


def video(forma, P, contraste=None, seed=0, radio=R):
    """P [n, nb, 2] (x, y) px verdaderos (NaN = no está). contraste [n, nb] (default CONTRASTE).

    Devuelve [n, alto, ancho] float32 en cuentas crudas (sin restar fondo).
    """
    rng = np.random.default_rng(seed)
    n, nb, _ = P.shape
    F = sint.dibujar(forma, np.zeros((0, 2)), radio)
    c = np.full((n, nb), CONTRASTE) if contraste is None else np.asarray(contraste, float)
    out = np.empty((n,) + tuple(forma), np.float32)
    for i in range(n):
        ideal = F.copy()
        for cv in np.unique(c[i]):
            sel = (c[i] == cv) & np.isfinite(P[i, :, 0])
            if sel.any() and cv > 0:
                ideal += sint.dibujar(forma, P[i, sel], radio, contraste=cv, fondo=F) - F
        out[i] = sint.camara(ideal, NIVEL, rng=rng)
    return out


def etiquetas(r, P, tol=1.0):
    """Bead verdadera de cada fila medida de r (−1 si ninguna a menos de tol px)."""
    lab = np.full(r['cuadro'].size, -1)
    for k, (i, x, y) in enumerate(zip(r['cuadro'], r['x'], r['y'])):
        if not np.isfinite(x):
            continue
        d = np.hypot(P[i, :, 0] - x, P[i, :, 1] - y)
        if np.nanmin(d) < tol:
            lab[k] = int(np.nanargmin(d))
    return lab


def pureza(r, P, tol=1.0):
    """{particula: (pureza, bead mayoritaria, filas medidas)} sobre las filas medidas."""
    lab = etiquetas(r, P, tol)
    out = {}
    for p in np.unique(r['particula']):
        s = (r['particula'] == p) & ~r['interpolado'] & np.isfinite(r['x'])
        v = lab[s]
        maj = np.bincount(v[v >= 0]).argmax() if (v >= 0).any() else -1
        out[int(p)] = ((v == maj).mean() if v.size else np.nan, int(maj), int(s.sum()))
    return out


def browniano(rng, n, x0, sigma, deriva=(0.0, 0.0)):
    """[n, nb, 2] paseo al azar (sigma px por paso y por eje) + deriva común (px/cuadro)."""
    x0 = np.asarray(x0, float)
    pasos = rng.normal(0, sigma, (n, len(x0), 2))
    pasos[0] = 0
    return x0[None] + np.cumsum(pasos, axis=0) + np.arange(n)[:, None, None] * np.asarray(deriva)


# ---------------------------------------------------------------- detección
@pytest.mark.parametrize('prep', ['crudo', 'dividido', 'restado'])
def test_detectar_todas_precision_y_normalizacion(prep):
    """Mismas detecciones con el cuadro crudo, dividido por el fondo (~1) o con fondo restado (~0)."""
    forma = (160, 200)
    P = np.array([[[40.3, 35.7], [120.6, 50.2], [170.1, 120.9], [60.8, 125.4]]] * 3)
    V = video(forma, P, seed=1)
    F = NIVEL * sint.dibujar(forma, np.zeros((0, 2)), R)
    V = {'crudo': V, 'dividido': V / F, 'restado': V - F}[prep]
    det = en.detectar_todas(lambda i: V[i], 3, R)
    assert det['cuadro'].size == 12 and list(det['n_por_cuadro']) == [4, 4, 4]
    d = np.hypot(det['x'][:, None] - P[0, None, :, 0], det['y'][:, None] - P[0, None, :, 1]).min(axis=1)
    assert d.max() < 0.3 and np.sqrt((d ** 2).mean()) < 0.15
    assert det['z'].min() > 15


def test_detectar_todas_cuadro_vacio_no_inventa():
    forma = (128, 160)
    V = video(forma, np.full((3, 1, 2), np.nan), seed=2)
    det = en.detectar_todas(lambda i: V[i], 3, R)
    assert det['cuadro'].size == 0


# ---------------------------------------------------------------- enlace (arrays)
def test_enlazar_memoria_y_busqueda_que_crece_con_el_hueco():
    """Hueco de 3 cuadros (= memoria) no corta; de 4 sí. La búsqueda crece como sqrt(hueco + 1)."""
    s = 2.0
    c = np.r_[np.arange(0, 10), np.arange(13, 20), np.arange(24, 30)]
    x = np.zeros(c.size)
    x[c >= 13] = 1.9 * s            # tras 3 cuadros sin verla: permitido hasta 2 s
    x[c >= 24] += 0.5
    det = dict(cuadro=c, x=10 + x, y=np.full(c.size, 10.0))
    p = en.enlazar(det, None, s, 3, prediccion_cuadros=0)
    assert set(p[c < 20]) == {0} and set(p[c >= 24]) == {1}
    # el mismo salto con hueco de 2 cuadros (búsqueda sqrt(3) s = 1.73 s < 1.9 s): se corta
    c2 = np.r_[np.arange(0, 10), np.arange(12, 20)]
    det2 = dict(cuadro=c2, x=10 + np.where(c2 >= 12, 1.9 * s, 0.0), y=np.full(c2.size, 10.0))
    p2 = en.enlazar(det2, None, s, 3, prediccion_cuadros=0)
    assert len(set(p2)) == 2


def test_enlazar_tiempo_real_cuadro_perdido_agranda_la_busqueda():
    """Un cuadro perdido de la cámara (salto en t) cuenta como Δt = 2 dt aunque los índices sean seguidos."""
    c = np.arange(6)
    t = np.r_[0, 1, 2, 4, 5, 6] * 0.01
    x = np.array([0, 0, 0, 1.3, 1.3, 1.3]) + 10       # 1.3 > s = 1, pero < sqrt(2) s
    det = dict(cuadro=c, x=x, y=np.full(6, 10.0))
    assert len(set(en.enlazar(det, t, 1.0, 0, 0))) == 1
    assert len(set(en.enlazar(det, None, 1.0, 0, 0))) == 2


def test_enlazar_asignacion_optima_no_codiciosa():
    """Lo más cercano primero dejaría una sin enlazar; la asignación óptima enlaza las dos."""
    det = dict(cuadro=np.array([0, 0, 1, 1]), x=np.array([0.0, 3.0, 2.0, 5.5]) + 10,
               y=np.full(4, 10.0))
    p = en.enlazar(det, None, 3.0, 0, 0)
    assert p[2] == p[0] and p[3] == p[1] and p[0] != p[1]


# ---------------------------------------------------------------- postprocesar (arrays)
def test_interpolacion_exacta_lineal_en_t_real_y_marcada():
    n = 50
    dt = 0.01
    t = np.arange(n) * dt
    t[11:] += dt                       # cuadro perdido de la cámara entre 10 y 11
    xv, yv = 3 + 70 * t, 5 - 20 * t   # px
    c = np.setdiff1d(np.arange(n), np.r_[10, 11, 30, 31, 32, 33, 34])
    crit = dict(interpolar_max_cuadros=3, memoria_cuadros=6, largo_min_cuadros=5)
    r = en.postprocesar(c, np.zeros(c.size, int), xv[c], yv[c], np.ones(c.size), t, R, crit)
    assert list(r['cuadro']) == list(np.setdiff1d(np.arange(n), np.arange(30, 35)))
    assert set(r['particula']) == {0}
    ii = r['interpolado']
    assert list(r['cuadro'][ii]) == [10, 11]
    np.testing.assert_allclose(r['x'][ii], xv[[10, 11]], atol=1e-12)
    np.testing.assert_allclose(r['y'][ii], yv[[10, 11]], atol=1e-12)
    assert np.isnan(r['z'][ii]).all() and np.isfinite(r['z'][~ii]).all()
    np.testing.assert_array_equal(r['x'][~ii], xv[r['cuadro'][~ii]])
    assert r['hueco_previo'][r['cuadro'] == 35][0] == 5 and r['hueco_previo'].sum() == 5
    res = r['resumen']
    assert res['n_interpolados'][0] == 2 and res['n_cuadros'][0] == n - 5
    assert res['n_medidos'][0] == n - 7


def test_hueco_mayor_que_memoria_corta():
    c = np.r_[np.arange(0, 30), np.arange(33, 60), np.arange(65, 90)]     # huecos de 3 y 5
    r = en.postprocesar(c, np.zeros(c.size, int), c * 0.1, c * 0.0, None, None, R,
                        dict(memoria_cuadros=3, largo_min_cuadros=5))
    assert list(r['resumen']['primer_cuadro']) == [0, 65]
    assert list(r['resumen']['ultimo_cuadro']) == [59, 89]
    assert r['resumen']['cortes'] == 1 and r['resumen']['cortes_hueco'] == 1


def cruce(n=81, v=1.5, ocultar=()):
    """Dos partículas en recta que se cruzan en (100, 100) en el cuadro n // 2 (arrays exactos)."""
    i = np.arange(n)
    s = (i - n // 2) * v
    A = np.column_stack([100 + s, np.full(n, 100.0)])
    B = np.column_stack([np.full(n, 100.0), 100 + s])
    d = np.hypot(*(A - B).T)
    mB = ~np.isin(i, ocultar)
    c = np.r_[i, i[mB]]
    p = np.r_[np.zeros(n, int), np.ones(mB.sum(), int)]
    xy = np.r_[A, B[mB]]
    return c, p, xy[:, 0], xy[:, 1], d


@pytest.mark.parametrize('accion', en.ACCIONES)
def test_colisiones_en_los_cuadros_correctos(accion):
    c, p, x, y, d = cruce()
    dcol = 2.5 * R
    esperado = np.nonzero(d < dcol)[0]
    crit = dict(accion_colision=accion, largo_min_cuadros=5)
    r = en.postprocesar(c, p, x, y, None, None, R, crit)
    assert r['resumen']['colisiones'] == 2 * esperado.size
    if accion == 'marcar':
        assert r['cuadro'].size == c.size
        for q in (0, 1):
            s = r['particula'] == q
            assert list(r['cuadro'][s & r['colision']]) == list(esperado)
        assert np.isfinite(r['x']).all()
    elif accion == 'descartar':
        assert r['cuadro'].size == c.size
        assert np.isnan(r['x'][r['colision']]).all() and np.isfinite(r['x'][~r['colision']]).all()
        for q in (0, 1):
            assert list(r['cuadro'][(r['particula'] == q) & r['colision']]) == list(esperado)
        assert r['resumen']['n_medidos'].tolist() == [c.size // 2 - esperado.size] * 2
    else:
        assert not np.isin(r['cuadro'], esperado).any() and not r['colision'].any()
        assert list(r['particula_origen']) == sorted(r['particula_origen'])
        assert list(r['resumen']['particula_origen']) == [0, 0, 1, 1]
        assert list(r['resumen']['ultimo_cuadro']) == [esperado[0] - 1, 80, esperado[0] - 1, 80]
        assert list(r['resumen']['primer_cuadro']) == [0, esperado[-1] + 1] * 2
        assert r['resumen']['cortes_colision'] == 2


def test_colision_con_una_oculta_se_detecta_por_interpolacion():
    """Al cruzarse se funden en una detección: B falta 3 cuadros, A igual se marca ahí."""
    oc = (39, 40, 41)
    c, p, x, y, d = cruce(ocultar=oc)
    esperado = np.nonzero(d < 2.5 * R)[0]
    r = en.postprocesar(c, p, x, y, None, None, R, dict(largo_min_cuadros=5))
    sA = r['particula'] == 0
    assert list(r['cuadro'][sA & r['colision']]) == list(esperado)
    sB = r['particula'] == 1
    assert list(r['cuadro'][sB & r['colision']]) == [k for k in esperado if k not in oc]
    # con interpolación: B se rellena en el hueco y queda marcada
    r = en.postprocesar(c, p, x, y, None, None, R, dict(largo_min_cuadros=5, interpolar_max_cuadros=3))
    sB = r['particula'] == 1
    assert list(r['cuadro'][sB & r['colision']]) == list(esperado)
    assert list(r['cuadro'][sB & r['interpolado']]) == list(oc)
    # fusión más larga que la memoria: B se corta en dos, A igual queda marcada en toda la fusión
    oc9 = tuple(range(36, 45))
    c9, p9, x9, y9, _ = cruce(ocultar=oc9)
    r9 = en.postprocesar(c9, p9, x9, y9, None, None, R, dict(largo_min_cuadros=5, memoria_cuadros=3))
    assert r9['resumen']['cortes_hueco'] == 1
    sA = r9['particula'] == 0
    assert list(r9['cuadro'][sA & r9['colision']]) == list(esperado)
    # y con 'cortar' A también se parte ahí (nada fundido queda en una trayectoria)
    r9 = en.postprocesar(c9, p9, x9, y9, None, None, R,
                         dict(largo_min_cuadros=5, memoria_cuadros=3, accion_colision='cortar'))
    assert not np.isin(r9['cuadro'], esperado).any()
    # descartar: no se interpola a través de la colisión
    r = en.postprocesar(c, p, x, y, None, None, R,
                        dict(largo_min_cuadros=5, interpolar_max_cuadros=3, accion_colision='descartar'))
    assert not r['interpolado'].any()


def test_cortas_se_descartan():
    largos = [5, 19, 20, 50]
    c = np.concatenate([np.arange(L) for L in largos])
    p = np.concatenate([np.full(L, k) for k, L in enumerate(largos)])
    x = 30.0 * p + 0.01 * c
    r = en.postprocesar(c, p, x, np.zeros(c.size), None, None, R, dict(largo_min_cuadros=20))
    assert sorted(set(r['particula'])) == [0, 1]
    assert list(r['resumen']['particula_origen']) == [2, 3]
    assert list(r['resumen']['n_cuadros']) == [20, 50]
    assert r['resumen']['descartadas_cortas'] == 2


# ---------------------------------------------------------------- de punta a punta (video)
def test_libres_entran_y_salen_con_deriva_pureza():
    """8 beads brownianas (0.5 px/paso) con deriva 0.8 px/cuadro: entran por la izquierda y salen
    por la derecha. Más una espuria de 10 cuadros. Cada bead = UNA trayectoria pura."""
    rng = np.random.default_rng(3)
    n, forma = 120, (190, 256)
    x0 = [(xx, yy) for xx in (-40, 40, 120, 200) for yy in (50, 140)]
    P = browniano(rng, n, x0, 0.5, deriva=(0.8, 0.0))
    esp = np.full((n, 1, 2), np.nan)
    esp[50:60, 0] = (160.0, 95.0)
    PP = np.concatenate([P, esp], axis=1)
    V = video(forma, PP, seed=4)
    fr = []
    r = en.seguir_todas(lambda i: V[i], n, np.arange(n) * 0.01, R, {}, progreso=fr.append)
    assert fr[-1] == 1.0 and np.all(np.diff(fr) >= 0)
    pz = pureza(r, PP)
    assert len(pz) == 8
    assert all(v[0] == 1.0 for v in pz.values())
    assert sorted(v[1] for v in pz.values()) == list(range(8))          # una por bead, sin la espuria
    assert r['resumen']['descartadas_cortas'] >= 1
    lab = etiquetas(r, PP)
    err = np.hypot(r['x'] - PP[r['cuadro'], lab, 0], r['y'] - PP[r['cuadro'], lab, 1])
    assert np.sqrt(np.mean(err ** 2)) < 0.15
    # cada bead se ve en todos los cuadros en que está entera adentro (sin huecos)
    h = tr.lado_recorte_sugerido(R) // 2
    for p, (_, b, m) in pz.items():
        dentro = ((P[:, b, 0] >= h + 1) & (P[:, b, 0] <= forma[1] - h - 2)).sum()
        assert m >= dentro - 2
    assert not r['colision'].any() and r['resumen']['cortes'] == 0


def _parpadeo(n=100):
    rng = np.random.default_rng(5)
    P = browniano(rng, n, [(60, 60), (150, 70), (100, 140)], 0.3)
    c = np.full((n, 3), CONTRASTE)
    c[30:32, 0] = 0.01                   # 2 cuadros casi invisible
    c[60:66, 0] = 0.01                   # 6 cuadros
    return P, video((200, 220), P, contraste=c, seed=6)


def test_parpadeo_memoria_corta_o_no():
    n = 100
    P, V = _parpadeo(n)
    t = np.arange(n) * 0.01
    det = en.detectar_todas(lambda i: V[i], n, R)
    assert list(det['n_por_cuadro'][[30, 31, 60, 65]]) == [2, 2, 2, 2]     # se dejó de ver
    r = en.seguir_todas(lambda i: V[i], n, t, R, dict(memoria_cuadros=3))
    pz = pureza(r, P)
    assert len(pz) == 4 and all(v[0] == 1.0 for v in pz.values())
    s0 = [p for p, v in pz.items() if v[1] == 0]
    assert len(s0) == 2                                                   # hueco de 6 > 3 cortó
    a, b = sorted(s0, key=lambda p: r['cuadro'][r['particula'] == p][0])
    ca, cb = r['cuadro'][r['particula'] == a], r['cuadro'][r['particula'] == b]
    assert ca[0] == 0 and ca[-1] == 59 and cb[0] == 66 and cb[-1] == n - 1
    assert r['hueco_previo'][(r['particula'] == a) & (r['cuadro'] == 32)][0] == 2   # el de 2 no cortó
    assert r['resumen']['cortes'] == 0           # el corte lo hizo enlazar (memoria), no un hueco
    # memoria 6 → una sola trayectoria; con interpolar 2 → el hueco de 2 se rellena y se marca
    r = en.seguir_todas(lambda i: V[i], n, t, R, dict(memoria_cuadros=6, interpolar_max_cuadros=2))
    pz = pureza(r, P)
    assert len(pz) == 3 and all(v[0] == 1.0 for v in pz.values())
    p0 = [p for p, v in pz.items() if v[1] == 0][0]
    s = r['particula'] == p0
    assert list(r['cuadro'][s & r['interpolado']]) == [30, 31]
    assert r['hueco_previo'][s & (r['cuadro'] == 66)][0] == 6
    err = np.hypot(r['x'][s & r['interpolado']] - P[[30, 31], 0, 0],
                   r['y'][s & r['interpolado']] - P[[30, 31], 0, 1])
    assert err.max() < 1.5


@pytest.mark.parametrize('accion', en.ACCIONES)
def test_cruce_en_video(accion):
    """Dos beads en recta (1.5 px/cuadro) que se cruzan en el cuadro 40 y se funden en UNA
    detección durante ~9 cuadros (d < 1.5 R), más que la memoria (3).

    Lo que se garantiza: las colisiones se marcan en los cuadros correctos (truth d < 2.5 R, con
    ±3 px de tolerancia por el sesgo de la posición fundida) y, con 'cortar', todas las
    trayectorias son puras. Con 'marcar'/'descartar' los ids ANTES del cruce son puros; después
    de una fusión más larga que la memoria el id puede pasar a la otra bead (por eso existe
    'cortar'; las filas dudosas quedan con colision = 1).
    """
    n = 81
    c, p, x, y, d = cruce(n)
    P = np.stack([np.column_stack([x[p == 0], y[p == 0]]), np.column_stack([x[p == 1], y[p == 1]])], axis=1)
    V = video((200, 200), P, seed=7)
    dcol = 2.5 * R
    r = en.seguir_todas(lambda i: V[i], n, None, R, dict(accion_colision=accion))
    pz = pureza(r, P, tol=1.5)
    lab = etiquetas(r, P, tol=1.5)
    lejos = d[r['cuadro']] > dcol + 3
    assert (lab[lejos] >= 0).all()                     # lejos de la colisión: posiciones buenas
    if accion != 'cortar':
        antes = lejos & (r['cuadro'] < n // 2)
        assert sorted(set(lab[antes])) == [0, 1]
        for q in set(r['particula'][antes]):
            assert len(set(lab[antes & (r['particula'] == q)])) == 1
        fl = r['cuadro'][r['colision']]
        assert (d[fl] < dcol + 3).all()                # nada marcado lejos
        visibles = np.nonzero((d < dcol - 3) & (d > 2 * R))[0]   # cerca pero todavía separadas
        for i in visibles:
            assert (r['colision'] & (r['cuadro'] == i)).sum() == 2
        fundidas = np.nonzero(d < 1.0 * R)[0]          # una sola detección: igual marcada
        assert all((r['colision'] & (r['cuadro'] == i)).any() for i in fundidas)
        if accion == 'descartar':
            assert np.isnan(r['x'][r['colision']]).all() and np.isfinite(r['x'][~r['colision']]).all()
    else:
        assert len(pz) >= 4 and all(v[0] == 1.0 for v in pz.values())
        assert (d[r['cuadro']] >= dcol - 3).all()
        assert sorted({v[1] for v in pz.values()}) == [0, 1]


def test_stokes_prediccion_mantiene_id_en_un_hueco():
    """Beads pegadas que cruzan el campo a 3 px/cuadro; una no se ve 3 cuadros (memoria = 3).
    Con predicción sigue con el mismo id; sin predicción (búsqueda sqrt(4)·5 = 10 px < 12 px) se corta."""
    n, forma = 70, (160, 240)
    v = 3.0
    x0 = np.array([(-30.0, 40.0), (40.0, 40.0), (110.0, 40.0), (10.0, 110.0), (80.0, 110.0)])
    P = x0[None] + np.arange(n)[:, None, None] * np.array([v, 0.0])
    c = np.full((n, len(x0)), CONTRASTE)
    c[30:33, 1] = 0.01
    V = video(forma, P, contraste=c, seed=8)
    t = np.arange(n) * 0.01
    for pred, n_tray_b1 in ((10, 1), (0, 2)):
        r = en.seguir_todas(lambda i: V[i], n, t, R, dict(prediccion_cuadros=pred, largo_min_cuadros=10))
        pz = pureza(r, P)
        assert all(q[0] == 1.0 for q in pz.values())
        assert sum(q[1] == 1 for q in pz.values()) == n_tray_b1
        assert {q[1] for q in pz.values()} == set(range(len(x0)))


# ---------------------------------------------------------------- desde_seguir
def test_desde_seguir_reproduce_seguir_sin_huecos():
    rng = np.random.default_rng(9)
    n, forma = 60, (120, 160)
    centros = np.array([(50.3, 60.7), (110.6, 55.2)])
    P = centros[None] + rng.normal(0, 0.4, (n, 2, 2))
    V = video(forma, P, seed=10)
    s = tr.seguir(lambda i: V[i].astype(float), n, centros, radio_px=R)
    assert s['valido'].all()
    t = np.arange(n) * 0.005
    r = en.desde_seguir(s, t, R, {})
    assert sorted(set(r['particula'])) == [0, 1]
    for b in (0, 1):
        sel = r['particula'] == b
        np.testing.assert_array_equal(r['cuadro'][sel], np.arange(n))
        np.testing.assert_array_equal(r['x'][sel], s['x'][:, b])
        np.testing.assert_array_equal(r['y'][sel], s['y'][:, b])
        np.testing.assert_array_equal(r['z'][sel], s['contraste_z'][:, b])
    assert not r['interpolado'].any() and not r['colision'].any() and not r['hueco_previo'].any()
    assert list(r['resumen']['n_cuadros']) == [n, n] and r['resumen']['cortes'] == 0


def test_desde_seguir_huecos_memoria_e_interpolacion():
    n = 80
    x = np.linspace(20, 30, n)[:, None]
    y = np.full((n, 1), 40.0)
    x[20:22] = np.nan                    # 2 cuadros: no corta
    x[50:55] = np.nan                    # 5 cuadros > memoria 3: corta
    s = dict(x=x, y=np.where(np.isnan(x), np.nan, y))
    r = en.desde_seguir(s, None, R, dict(memoria_cuadros=3, interpolar_max_cuadros=2))
    assert list(r['resumen']['primer_cuadro']) == [0, 55]
    assert list(r['resumen']['ultimo_cuadro']) == [49, n - 1]
    assert list(r['cuadro'][r['interpolado']]) == [20, 21]
    np.testing.assert_allclose(r['x'][r['interpolado']], np.linspace(20, 30, n)[[20, 21]], atol=1e-12)
    assert r['resumen']['cortes_hueco'] == 1
