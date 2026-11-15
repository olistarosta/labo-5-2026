"""Tiempo: t_i con error < 0.1 ms y TODOS los cuadros perdidos detectados (criterio del encargo).

Condiciones: jitter del host exponencial de 2 ms + colas de 10-50 ms (FIFO), reloj de cámara
13.7 ppm más rápido que el nominal, cuadros perdidos al azar. Los números medidos se imprimen
(pytest -s) para el informe.
"""
import numpy as np
import pytest

from pinzas import io, sintetico, tiempo

CASOS = [(3000, 100.0, 2e-3), (20000, 200.0, 5e-4), (2000, 35.0, 20e-3)]   # n, fps, W


def _sin_hw(tm):
    tn = dict(tm)
    tn['t_camara'] = np.full(len(tm['t_camara']), np.nan)
    return tn


@pytest.mark.parametrize('n,fps,W', CASOS)
@pytest.mark.parametrize('seed', [0, 1, 2])
def test_con_timestamp_hardware(n, fps, W, seed):
    tm = sintetico.tiempos_camara(n, fps, W, prob_perdido=0.003, rng=seed)
    r = tiempo.reconstruir(tm)                               # tick_hz ajustado contra el host
    err = np.max(np.abs(r['t'] - tm['t_verdad']))
    print(f'\n[hw] n={n} fps={fps} seed={seed}: error max {err * 1e6:.2f} us, '
          f'perdidos {r["cuadros_perdidos"]}/{len(tm["contadores_perdidos"])}, '
          f'tick ajustado {r["informe"]["tick_hz_ajustado"]:.1f} Hz')
    assert err < 1e-4
    # todos los perdidos, y en el lugar justo
    faltan = np.concatenate([np.arange(a, b + 1) for a, b in r['informe']['posiciones_perdidos']]
                            or [np.zeros(0, int)])
    np.testing.assert_array_equal(faltan, tm['contadores_perdidos'])
    assert r['cuadros_perdidos'] == len(tm['contadores_perdidos'])
    assert r['t'][0] == 0.0
    # el reloj ajustado ve los 13.7 ppm del reloj real (medido con el host)
    assert abs(r['informe']['tick_hz_ajustado'] / 1.0000137e9 - 1) < 1e-6


def test_con_tick_nominal_el_error_es_el_de_los_ppm():
    """Si se usa el tick_hz nominal (lo que diría el SDK), el error crece como ppm × duración.

    Esto NO cumple 0.1 ms en grabaciones largas: queda documentado y el informe lo muestra
    (ppm_reloj, diferencia_max_vs_ajuste_s).
    """
    tm = sintetico.tiempos_camara(3000, 100.0, 2e-3, rng=5)
    r = tiempo.reconstruir(tm, tick_hz=1e9)
    err = np.max(np.abs(r['t'] - tm['t_verdad']))
    T = tm['t_verdad'][-1]
    print(f'\n[nominal] error max {err * 1e6:.1f} us en {T:.1f} s; ppm_reloj '
          f'{r["informe"]["ppm_reloj"]:.2f}; diferencia vs ajuste '
          f'{r["informe"]["diferencia_max_vs_ajuste_s"] * 1e6:.1f} us')
    assert err == pytest.approx(13.7e-6 / 1.0000137 * T, rel=0.02)
    assert r['informe']['ppm_reloj'] == pytest.approx(13.7, abs=0.2)
    assert r['informe']['diferencia_max_vs_ajuste_s'] == pytest.approx(err, rel=0.05)
    assert r['informe']['tick_hz_usado'] == 1e9


@pytest.mark.parametrize('n,fps,W', CASOS)
def test_sin_timestamp_hardware(n, fps, W):
    """Método de fluidos (recta local ±15). Criterio ~1 ms: se mide y se reporta.

    t_i = ŷ_i - ŷ_0, así que el error del primer ajuste se arrastra a TODA la grabación como
    una constante; por eso se mide también el error salvo una constante (lo que importa para
    los retardos). Cerca de una cola del host el error local llega a varios ms.
    """
    tm = sintetico.tiempos_camara(n, fps, W, prob_perdido=0.003, rng=11)
    r = tiempo.reconstruir(_sin_hw(tm))
    d = r['t'] - tm['t_verdad']
    d0 = d - np.median(d)
    rg = tiempo.reconstruir(_sin_hw(tm), ajuste_sin_hw='global')
    eg = np.abs(rg['t'] - tm['t_verdad'])
    print(f'\n[host local] n={n} fps={fps}: rms {np.sqrt(np.mean(d ** 2)) * 1e3:.3f} ms, '
          f'max {np.abs(d).max() * 1e3:.2f} ms; salvo constante: rms '
          f'{np.sqrt(np.mean(d0 ** 2)) * 1e3:.3f} ms, p95 {np.percentile(np.abs(d0), 95) * 1e3:.3f} ms, '
          f'max {np.abs(d0).max() * 1e3:.2f} ms | [host global] max {eg.max() * 1e6:.1f} us')
    assert r['informe']['metodo'] == 'host'
    assert r['cuadros_perdidos'] == len(tm['contadores_perdidos'])
    assert np.sqrt(np.mean(d0 ** 2)) < 1.5e-3
    assert np.percentile(np.abs(d0), 95) < 2.5e-3
    assert eg.max() < 1e-4


def test_contador_que_se_repite_o_retrocede_se_reporta():
    tm = sintetico.tiempos_camara(500, 100.0, 2e-3, prob_perdido=0.0, rng=3)
    c = tm['contador_camara'].copy()
    c[100] = c[99]          # leído antes de actualizarse: saltos 0 y 2
    c[300:] -= 5            # retrocede
    tm['contador_camara'] = c
    r = tiempo.reconstruir(tm)
    inf = r['informe']
    assert inf['patron_0_2'] == 1
    assert [a[0] for a in inf['anomalias_contador']] == [99, 299]
    # el 0-2 aparece como un hueco de 1 (no se repara: se reporta)
    assert list(r['huecos']) == [100] and r['cuadros_perdidos'] == 1
    # los tramos cortan en huecos y anomalías, y cubren todo sin superponerse
    tr = r['tramos']
    assert tr[0, 0] == 0 and tr[-1, 1] == 499
    assert np.all(tr[1:, 0] == tr[:-1, 1] + 1)
    assert set(tr[:-1, 1]) == {99, 100, 299}
    # el tiempo de hardware no depende del contador
    assert np.max(np.abs(r['t'] - tm['t_verdad'])) < 1e-4


def test_tramos_y_estadistica():
    tm = sintetico.tiempos_camara(4000, 100.0, 2e-3, prob_perdido=0.004, rng=8)
    r = tiempo.reconstruir(tm)
    inf = r['informe']
    assert len(r['tramos']) == len(r['huecos']) + 1
    assert inf['fps_medido'] == pytest.approx(1 / tm['periodo_real_s'], rel=1e-6)
    assert inf['dt_medio_s'] == pytest.approx(tm['periodo_real_s'], rel=1e-5)
    assert inf['dt_desvio_s'] < 5e-6                         # jitter de cámara 1 us
    assert inf['dt_max_s'] >= 2 * tm['periodo_real_s'] * 0.999  # hay huecos
    assert inf['residuo_max_host_s'] > 0.01                  # las colas se ven en el residuo
    assert inf['residuo_mediana_host_s'] == pytest.approx(2e-3 * np.log(2), rel=0.2)
    assert len(r['dt']) == len(r['t']) - 1


def test_t_camara_parcialmente_vacio_queda_nan():
    tm = sintetico.tiempos_camara(300, 100.0, 2e-3, prob_perdido=0.0, rng=2)
    tc = tm['t_camara'].astype(float)
    tc[[10, 11, 200]] = np.nan
    tm['t_camara'] = tc
    r = tiempo.reconstruir(tm)
    assert np.isnan(r['t'][[10, 11, 200]]).all()
    assert r['informe']['t_camara_faltantes'] == 3
    ok = np.isfinite(r['t'])
    assert np.max(np.abs(r['t'][ok] - tm['t_verdad'][ok])) < 1e-4


def test_informe_es_serializable_a_json(tmp_path):
    tm = sintetico.tiempos_camara(500, 100.0, 2e-3, rng=1)
    r = tiempo.reconstruir(tm)
    io.escribir_json(str(tmp_path / 'informe.json'), r['informe'])
    assert io.leer_json(str(tmp_path / 'informe.json'))['cuadros_perdidos'] == r['cuadros_perdidos']


def test_reconstruir_desde_disco(tmp_path):
    """Ida y vuelta por tiempos.csv (t_camara int64 → float64 sin perder ticks)."""
    car = sintetico.grabar(str(tmp_path), 'tiempo', 'libre', n=400, fps=100.0, W=2e-3,
                           beads=[], forma=(16, 16), prob_perdido=0.01, seed=4)
    g = io.leer_grabacion(car)
    r = tiempo.reconstruir(g.tiempos)
    t_v = np.asarray(g.meta['verdad']['t_inicio_s'])
    assert np.max(np.abs(r['t'] - t_v)) < 1e-4
    assert r['cuadros_perdidos'] == len(g.meta['verdad']['contadores_perdidos'])


def test_estadistica_en_vivo():
    tm = sintetico.tiempos_camara(1000, 100.0, 2e-3, prob_perdido=0.01, rng=6)
    e = tiempo.estadistica_en_vivo(tm['contador_camara'], tm['t_camara'], tm['t_host_ns'],
                                   tick_hz=1e9, n_ultimos=200)
    assert e['fuente'] == 'camara'
    assert e['fps'] == pytest.approx(100.0, rel=1e-4)        # nominal: 13.7 ppm de diferencia
    assert e['fps_host'] == pytest.approx(100.0, rel=0.03)    # jitter + colas del host
    assert e['perdidos'] == len(tm['contadores_perdidos'])
    e2 = tiempo.estadistica_en_vivo(tm['contador_camara'], None, tm['t_host_ns'], n_ultimos=200)
    assert e2['fuente'] == 'host' and e2['fps'] == e['fps_host']
    assert tiempo.estadistica_en_vivo([5])['n'] == 1


def test_envolvente_inferior_recupera_recta_bajo_la_nube():
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 1e10, 5000))
    y = 3.0 + 1.00002e-9 * x + rng.exponential(1e-3, x.size)
    y[::97] += 0.05                                          # colas
    a, b = tiempo.envolvente_inferior(x, y)
    assert b == pytest.approx(1.00002e-9, rel=1e-7)
    assert a == pytest.approx(3.0, abs=2e-5)
