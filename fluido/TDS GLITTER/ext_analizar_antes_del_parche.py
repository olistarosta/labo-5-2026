"""
Analizador — la lógica de /analizar. Lo numérico y lo visual está en los nodos; esta
clase solo decide QUÉ se calcula en cada cuadro de TD y guarda los resultados.

  Cargar()          al elegir grabación o cámara: carpetas, tiempos, dt
  Marcar(u, v)      click en la imagen: regla, centro, ROI, colores, ventana de ejemplo
  InicioCuadro()    onFrameStart: elige par y pasada y cocina la pasada en orden
  FinCuadro()       onFrameEnd: si terminó un par de ANALIZAR TODO, lo escribe
  Analizar()        recorre todos los pares de todas las cámaras
  Detener()         corta y guarda lo hecho
  GuardarEtapas()   imágenes del par que se ve (detalle y desplazamiento por iteración)

Salida: una carpeta por análisis, dentro de una carpeta por grabación:
  data glitter/<grabación>/<fecha y hora del análisis>/
      parametros.json          calibración, colores, seguimiento, dt, método
      camara<L>_campo.csv      un vector por fila, por par
      camara<L>_centro.csv     centro del vórtice por par (partícula seguida o fijo)
      imagenes/camara<L>_par<k>/   detalle_iterN y desplazamiento_iterN, sin texto (GUARDAR IMÁGENES)
ANALIZAR TODO guarda también las imágenes del par que se estaba mirando. GUARDAR IMÁGENES
agrega imágenes al último análisis de esa grabación (o crea la carpeta si no hay).

PASO TEMPORAL.
  Cámara: la hora de cada cuadro viene redondeada (~4 ms) pero la cámara tiene un ritmo
  fijo, así que la velocidad usa un paso fijo: dt = duración / intervalos, exacto sobre
  muchos cuadros. Un cuadro perdido se ve como un hueco y ese par usa 2*dt.
  Simulación: su reloj es exacto pero el ritmo de TD varía, así que cada par usa la
  diferencia de horas de sus dos cuadros (con dt fijo la velocidad salía hasta 20 % baja
  en los tramos en que TD iba más rápido).
"""

import datetime
import json
import os

import cv2
import numpy as np

CAMARAS = ('A', 'B')

# Orden de cocinado de una pasada. Todo lee 'memoria', que no cambia durante el
# cuadro, así que el resultado no depende de quién tire de la red ni en qué orden.
CADENA = ('a_denso', 'predictor_en_grilla', 'deformar_a', 'deformar_b',
          'correlacion', 'sumar', 'validar', 'suavizar', 'campo_suave')

# Vistas del proceso: guardan los datos de cada pasada en por_pasada (en este orden:
# la correlación usa el paso que guardó la grilla).
VISTAS_POR_PASADA = ('superposicion_ab', 'dibujar_grilla', 'dibujar_correlacion', 'dibujar_vectores')

LADO_IMAGEN = 2048      # px del lado de las imágenes exportadas (recortadas al ROI)
LADO_VENTANA = 1024     # px de las ventanas y el plano de correlación exportados


class Analizador:
    def __init__(self, ownerComp):
        self.ownerComp = ownerComp
        self.T = None           # hora de cada cuadro de la cámara activa, s
        self.Nint = None        # cuadros que abarca cada intervalo (1 salvo pérdidas)
        self.RelojExacto = False  # True en la simulación: dt de cada par = diferencia de horas
        self.Camaras = {}       # 'A'/'B' -> carpeta con cuadros/ y tiempos.csv
        self.Info = {}          # grabacion.json
        self.Cola = []          # (par, pasada) pendientes
        self.EnCurso = None     # (par, pasada) de este cuadro de TD
        self.Clave = None       # parámetros con que se calculó lo que se ve
        self.Corrida = None     # análisis completo en curso

    # ------------------------------------------------------------------ util
    def _npasadas(self):
        return int(self.ownerComp.par.Pasadas)

    def _separacion(self):
        """Cuántos cuadros hay entre los dos del par: se compara k con k + separación."""
        return int(self.ownerComp.par.Separacion)

    def _tiempo_par(self, k):
        """(dt del par, tiempo del centro del par) para el par k, k + separación.

        Cámara: dt fijo por cuadro por la cantidad de intervalos entre los dos cuadros (cuenta
        los cuadros perdidos). Simulación: diferencia exacta de horas.
        """
        k2 = k + self._separacion()
        if self.RelojExacto:
            return self.T[k2] - self.T[k], 0.5 * (self.T[k] + self.T[k2]) - self.T[0]
        dt = float(self.ownerComp.par.Dt)
        cum = np.r_[0, np.cumsum(self.Nint)]
        return (cum[k2] - cum[k]) * dt, 0.5 * (cum[k] + cum[k2]) * dt

    def _ventana(self, pasada):
        """Ventana de la pasada: 4v, 2v y después v (las pasadas extra repiten la final)."""
        return int(self.ownerComp.par.Ventana) * 2 ** max(0, 2 - pasada)

    def MenuPasadas(self):
        """Opciones de "Pasada que se ve" según la cantidad de pasadas."""
        p = self.ownerComp.par
        n = self._npasadas()
        p.Pasadaver.menuNames = [str(k + 1) for k in range(n)]
        p.Pasadaver.menuLabels = ['iteración %d (ventana %d px)' % (k + 1, self._ventana(k)) for k in range(n)]
        if int(p.Pasadaver.menuIndex) >= n or p.Pasadaver.eval() not in p.Pasadaver.menuNames:
            p.Pasadaver = str(n)

    def _say(self, msg):
        self.ownerComp.par.Estado = msg

    def _reproducir(self):
        """El cálculo avanza una pasada por cuadro de TD: con la línea de tiempo en pausa no hay
        onFrameStart y nada avanza, sin ningún error. Se pone en marcha al cargar y al analizar."""
        if not op('/').time.play:
            op('/').time.play = True

    def _grabacion(self):
        return self.ownerComp.par.Grabacion.eval().strip().rstrip('/\\')

    def _salida(self):
        """data glitter/<grabación>/: la carpeta con los análisis de esta grabación."""
        return os.path.normpath(os.path.join(project.folder, '..', 'data glitter',
                                             os.path.basename(self._grabacion())))

    def _nueva_carpeta_analisis(self):
        base = os.path.join(self._salida(), datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S'))
        carpeta, i = base, 2
        while os.path.exists(carpeta):
            carpeta, i = '%s_%d' % (base, i), i + 1
        os.makedirs(carpeta)
        return carpeta

    def _ultima_carpeta_analisis(self):
        raiz = self._salida()
        if not os.path.isdir(raiz):
            return None
        # Solo análisis terminados (con parametros.json): una carpeta cortada no cuenta.
        hechas = sorted(d for d in os.listdir(raiz) if os.path.isfile(os.path.join(raiz, d, 'parametros.json')))
        return os.path.join(raiz, hechas[-1]) if hechas else None

    # ----------------------------------------------------------------- carga
    def Cargar(self):
        comp = self.ownerComp
        self._reproducir()
        carpeta = self._grabacion()
        self.T = None
        self.Camaras = {}
        comp.par.Carpetacuadros = ''
        if not carpeta or not os.path.isdir(carpeta):
            comp.par.Info = ''
            self._say('1. Elegí en "Grabación" la carpeta que creó GRABAR.')
            return
        for L in CAMARAS:
            sub = os.path.join(carpeta, 'camara_' + L)
            if os.path.isdir(os.path.join(sub, 'cuadros')):
                self.Camaras[L] = sub
        if not self.Camaras and os.path.isdir(os.path.join(carpeta, 'cuadros')):
            self.Camaras['A'] = carpeta            # grabación vieja, de una sola cámara
        if not self.Camaras:
            comp.par.Info = ''
            self._say('Esa carpeta no tiene cuadros. Elegí la carpeta de una grabación.')
            return
        self.Info = {}
        jr = os.path.join(carpeta, 'grabacion.json')
        if os.path.exists(jr):
            with open(jr, encoding='utf-8') as f:
                self.Info = json.load(f)
        L = comp.par.Camara.eval()
        if L not in self.Camaras:
            comp.par.Camara = sorted(self.Camaras)[0]   # dispara otra carga, con esa cámara
            return
        sub = self.Camaras[L]
        t = np.loadtxt(os.path.join(sub, 'tiempos.csv'), delimiter=',', skiprows=1, ndmin=2)[:, 1]
        d = np.diff(t)
        d0 = np.median(d[d > 0]) if np.any(d > 0) else 1.0
        n = np.maximum(np.rint(d / d0), 1).astype(int)
        dt = (t[-1] - t[0]) / n.sum()
        n = np.maximum(np.rint(d / dt), 1).astype(int)      # afinar con el dt promediado
        dt = (t[-1] - t[0]) / n.sum()
        self.T, self.Nint = t, n
        cam = self.Info.get('camaras', {}).get(L, self.Info)
        self.RelojExacto = cam.get('fuente') == 'simulacion'
        comp.par.Carpetacuadros = os.path.join(sub, 'cuadros').replace('\\', '/')
        sep = self._separacion()
        comp.par.Npares = max(0, len(t) - sep)                  # pares (k, k + separación)
        comp.par.Dt = dt
        ultimo = max(0, len(t) - sep - 1)
        comp.par.Cuadro.normMax = max(1, ultimo)
        comp.par.Cuadro.max = ultimo
        comp.par.Cuadro.clampMax = True
        comp.par.Cuadro = min(int(comp.par.Cuadro), ultimo)
        sim = self._verdad(L)
        if sim and float(comp.par[L + 'roiradio']) == 0:
            # Grabación de la simulación: centro y ROI conocidos.
            centro = sim.get('centro_medio_px', sim.get('centro_px'))
            comp.par[L + 'centrox'], comp.par[L + 'centroy'] = centro
            comp.par[L + 'roiradio'] = sim['radio_sembrado_px'] + sim.get('deriva_px', 0.0)
        comp.par.Info = 'cámaras: %s | %s: %d cuadros, %.1f fps, %d perdidos | pares k y k+%d' % (
            ' y '.join(sorted(self.Camaras)), L, len(t), 1.0 / dt, int((n - 1).sum()), sep)
        aviso = self._ajustar_marcas(L)
        self.Clave = None
        if self.Corrida is None:
            self._say(aviso or self._pendiente() or 'Listo. Mirá el par en "Mostrar" y apretá ANALIZAR TODO.')

    def _ajustar_marcas(self, L):
        """Las marcas (regla, centro, ROI) son PÍXELES DE LA IMAGEN, así que solo valen
        para la resolución con la que se marcó.

        TouchDesigner sin licencia comercial limita los TOP a 1280 px, así que una
        grabación de 1920 x 1080 se analiza a 1280 x 720: marcas hechas sobre la imagen
        grande dejan el ROI corrido y la escala en mm/px 1.5 veces mal, sin ningún error
        a la vista. Cada cámara guarda con qué ancho se marcó y acá se reescala si hace
        falta. Devuelve el aviso para el estado, o None si no hubo que tocar nada.
        """
        comp = self.ownerComp
        ancho = float(comp.op('cuadro_a').width)
        if ancho <= 0:
            return None
        marcadas = float(comp.par[L + 'resolucion'])
        if marcadas <= 0:
            # Marcas viejas, de antes de que se guardara el ancho: si alguna cae afuera
            # de la imagen, se marcaron sobre la resolución original de la grabación.
            nativa = (self.Info.get('camaras', {}).get(L, self.Info)
                      .get('resolucion_px', [0, 0])[0])
            lejos = max(float(comp.par[L + 'puntoax']), float(comp.par[L + 'puntobx']),
                        float(comp.par[L + 'centrox']) + float(comp.par[L + 'roiradio']))
            marcadas = float(nativa) if (nativa and lejos > ancho) else ancho
        k = ancho / marcadas
        if abs(k - 1.0) < 1e-6:
            comp.par[L + 'resolucion'] = ancho
            return None
        for nombre in ('puntoax', 'puntoay', 'puntobx', 'puntoby', 'centrox', 'centroy', 'roiradio'):
            comp.par[L + nombre] = float(comp.par[L + nombre]) * k
        comp.par[L + 'resolucion'] = ancho
        return ('Las marcas de la cámara %s estaban en %d px de ancho y la imagen es de %d: '
                'se reescalaron x%.3f (mirá que el ROI siga bien).' % (L, marcadas, ancho, k))

    def _verdad(self, L):
        cam = self.Info.get('camaras', {}).get(L)
        if cam:
            return cam.get('simulacion')
        return self.Info.get('simulacion') if L == 'A' else None

    def _pendiente(self):
        """Qué falta, dicho en el orden en que conviene hacerlo."""
        p = self.ownerComp.par
        faltan = []
        for L in sorted(self.Camaras):
            if float(p[L + 'escala']) == 0:
                faltan.append('%s: regla (Marcar = Punto A, B + distancia)' % L)
            if float(p[L + 'roiradio']) == 0:
                faltan.append('%s: centro y ROI' % L)
        if not faltan:
            return ''
        return 'Calibración pendiente (sin ella las velocidades salen en px/s):\n  ' + '\n  '.join(faltan)

    # ------------------------------------------------------ click en la imagen
    def Marcar(self, u, v):
        """Click sobre la imagen (u, v en 0..1 del panel) -> pixel de la imagen.

        La imagen se dibuja ajustada al panel conservando la proporción: hay que
        descontar la banda que sobra o la calibración sale corrida en silencio. El
        tamaño del panel sale de SUS PARÁMETROS, nunca de su geometría calculada
        (leerla desde un script puede congelar el layout de TD).
        """
        comp = self.ownerComp
        modo = comp.par.Marcar.eval()
        if modo == 'nada':
            return
        panel = comp.op('panel_imagen')
        pw, ph = float(panel.par.w.eval()), float(panel.par.h.eval())
        img = comp.op('cuadro_a')
        iw, ih = float(img.width), float(img.height)
        esc = min(pw / iw, ph / ih)
        x = (u * pw - (pw - iw * esc) * 0.5) / esc
        y = (v * ph - (ph - ih * esc) * 0.5) / esc
        if not (0 <= x < iw and 0 <= y < ih):
            self._say('el click cayó fuera de la imagen')
            return
        p, L = comp.par, comp.par.Camara.eval()
        p[L + 'resolucion'] = iw          # las marcas valen para ESTE ancho de imagen
        if modo == 'puntoa':
            p[L + 'puntoax'], p[L + 'puntoay'] = x, y
            p.Marcar = 'puntob'
            self._say('%s: punto A = (%.0f, %.0f) px. Ahora marcá el punto B.' % (L, x, y))
        elif modo == 'puntob':
            p[L + 'puntobx'], p[L + 'puntoby'] = x, y
            p.Marcar = 'nada'
            self._say('%s: punto B = (%.0f, %.0f) px. Cargá la distancia real A-B en la página Cámara %s.' % (L, x, y, L))
        elif modo == 'centro':
            p[L + 'centrox'], p[L + 'centroy'] = x, y
            p.Marcar = 'roi'
            self._say('%s: centro = (%.0f, %.0f) px. Ahora marcá el borde del ROI.' % (L, x, y))
        elif modo == 'roi':
            p[L + 'roiradio'] = np.hypot(x - float(p[L + 'centrox']), y - float(p[L + 'centroy']))
            p.Marcar = 'nada'
            self._say('%s: radio del ROI = %.0f px.\n%s' % (L, float(p[L + 'roiradio']),
                                                         self._pendiente() or 'Calibración completa.'))
        elif modo in ('colorfondo', 'colorglitter', 'colorparticula'):
            rgb = self._color_en(x, y, modo)
            for c, val in zip('rgb', rgb):
                p[L + modo + c] = float(val)
            siguiente = {'colorfondo': 'colorglitter', 'colorglitter': 'colorparticula',
                         'colorparticula': 'nada'}[modo]
            p.Marcar = siguiente
            nombre = {'colorfondo': 'del fondo', 'colorglitter': 'del glitter',
                      'colorparticula': 'de la partícula'}[modo]
            pedir = {'colorglitter': '\nAhora hacé click sobre un destello de glitter.',
                     'colorparticula': '\nAhora hacé click sobre la partícula central.',
                     'nada': ''}[siguiente]
            self._say('%s: color %s = (%.2f, %.2f, %.2f)%s\n%s' % (
                L, nombre, rgb[0], rgb[1], rgb[2], pedir,
                comp.op('contraste_color').module.resumen(comp)))
        elif modo == 'ejemplo':
            p.Ejemplox, p.Ejemploy = x, y
            p.Marcar = 'nada'
            p.Vista = 'correlacion'
            self._say('Ventana de ejemplo en (%.0f, %.0f) px. Cambiá "Pasada que se ve" para comparar.' % (x, y))

    def _color_en(self, x, y, modo):
        """Color RGB alrededor del click.

        Fondo: mediana de 15x15 px. Glitter y partícula: los destellos son de pocos px,
        así que un promedio los diluye en el fondo; se toma el promedio del 25 % de los
        px de 9x9 más alejados del color del fondo.
        """
        comp = self.ownerComp
        img = comp.op('cuadro_a').numpyArray(delayed=False)[..., :3]
        h, w = img.shape[:2]
        i, j = int(round(x)), int(round(y))
        r = 7 if modo == 'colorfondo' else 4
        zona = img[max(0, j - r):min(h, j + r + 1), max(0, i - r):min(w, i + r + 1)].reshape(-1, 3)
        if modo == 'colorfondo':
            return np.median(zona, axis=0)
        fondo = np.array([float(comp.par['Colorfondo' + c]) for c in 'rgb'])
        dist = np.linalg.norm(zona - fondo, axis=1)
        lejos = zona[dist >= np.percentile(dist, 75)]
        return lejos.mean(axis=0)

    # ------------------------------------------------------ secuenciador
    def _clave(self):
        """Todo lo que cambia el resultado del par que se ve. Si cambia, se recalcula."""
        p = self.ownerComp.par
        nombres = ('Carpetacuadros', 'Cuadro', 'Separacion', 'Ventana', 'Pasadas', 'Imagen', 'Picomin',
                   'Centrox', 'Centroy', 'Roiradio',
                   'Seguirparticula', 'Umbral', 'Areaminima',
                   'Colorfondor', 'Colorfondog', 'Colorfondob', 'Colorglitterr', 'Colorglitterg',
                   'Colorglitterb', 'Colorparticular', 'Colorparticulag', 'Colorparticulab')
        return tuple(p[n].eval() for n in nombres)

    def _encolar(self, k):
        """Las pasadas del par k, desde la primera."""
        self.Cola = [(k, pasada) for pasada in range(self._npasadas())]

    def InicioCuadro(self):
        """onFrameStart: elige par y pasada, guarda la memoria y cocina la pasada."""
        comp, bucle = self.ownerComp, self.ownerComp.op('bucle')
        if self.T is None:
            return
        if not self.Cola and self.Corrida is None:
            clave = self._clave()
            if clave != self.Clave:
                self.Clave = clave
                comp.op('por_pasada').module.limpiar()
                comp.par.Referenciax, comp.par.Referenciay = comp.par.Centrox.eval(), comp.par.Centroy.eval()
                self._encolar(int(comp.par.Cuadro))
        if not self.Cola:
            self.EnCurso = None
            return
        self.EnCurso = self.Cola.pop(0)
        k, pasada = self.EnCurso
        comp.par.Paractual, comp.par.Pasada = k, pasada
        # 0) Partícula central: una vez por par, con los dos cuadros ya cargados.
        if pasada == 0:
            comp.op('centro_particula').cook(force=True)
        # 1) MEMORIA: el feedback devuelve campo_suave tal como quedó al final del
        #    cuadro anterior, o sea la pasada anterior. Se guarda en memoria, que es
        #    lo único que leen los nodos de esta pasada.
        bucle.op('campo_previo').cook(force=True)
        memoria = bucle.op('memoria')
        memoria.par.activepulse.pulse()
        memoria.cook(force=True)
        # 2) La pasada. Se cocina ACÁ, al inicio del cuadro, y no en onFrameEnd:
        #    medido, TD toma el objetivo del feedback al final de la fase de
        #    cocinado, ANTES de onFrameEnd. Cada nodo una sola vez por cuadro.
        cuadro = absTime.frame
        for nombre in CADENA:
            nodo = bucle.op(nombre)
            if nodo.cookAbsFrame != cuadro:
                nodo.cook(force=True)
        # 3) Las vistas del proceso guardan la imagen de esta pasada (no tocan el bucle).
        #    En ANALIZAR TODO no se dibujan: cuestan más que el PIV entero y no se miran
        #    (la pantalla muestra el cuadro original, que además deja ver cómo avanza).
        #    Al terminar, el par que se esté mirando se recalcula con todas sus vistas.
        if self.Corrida is None:
            for nombre in VISTAS_POR_PASADA:
                comp.op(nombre).cook(force=True)

    def FinCuadro(self):
        """onFrameEnd: si era la última pasada de un par de ANALIZAR TODO, lo escribe."""
        if self.EnCurso is None:
            return
        k, pasada = self.EnCurso
        if pasada == self._npasadas() - 1 and self.Corrida is not None:
            self._escribir(k)
            self._siguiente()

    # --------------------------------------------------- análisis completo
    def Analizar(self):
        if self.Corrida is not None:
            return
        self._reproducir()
        if not self.Camaras:
            self.Cargar()
            if not self.Camaras:
                return
        carpeta = self._nueva_carpeta_analisis()
        # Las imágenes del par que se estaba mirando, si ya está calculado (se pierden al empezar).
        etapas = ''
        if (not self.Cola and self.EnCurso is None
                and self.ownerComp.op('por_pasada').module.leer('suave', self._npasadas() - 1) is not None):
            etapas = '\nimágenes del par %d en imagenes/' % int(self.ownerComp.par.Cuadro)
            self._exportar_etapas(carpeta)
        self.Corrida = {'base': os.path.basename(self._grabacion()), 'camaras': sorted(self.Camaras),
                        'hechas': {}, 'inicio': datetime.datetime.now(), 'csv': None,
                        'salida': carpeta, 'etapas': etapas}
        self.ownerComp.par.Corriendo = True     # mientras corre no se dibujan las vistas
        self._empezar_camara()

    def Detener(self):
        if self.Corrida is None:
            return
        self.Corrida['detenido'] = True
        self._cerrar_camara()
        self._cerrar()

    def _empezar_camara(self):
        comp, r = self.ownerComp, self.Corrida
        L = r['camaras'][len(r['hechas'])]
        if comp.par.Camara.eval() != L:
            comp.par.Camara = L
        self.Cargar()                          # tiempos y carpeta de esta cámara
        ruta = os.path.join(r['salida'], 'camara%s_campo.csv' % L)
        r['csv'] = open(ruta, 'w', encoding='utf-8', newline='')
        r['csv'].write('frame,t_s,x_px,y_px,x_mm,y_mm,u_mm_s,v_mm_s,valido\n')
        ruta_c = os.path.join(r['salida'], 'camara%s_centro.csv' % L)
        r['csv_centro'] = open(ruta_c, 'w', encoding='utf-8', newline='')
        r['csv_centro'].write('frame,t_s,xc_px,yc_px,xc_mm,yc_mm,seguido,area_px\n')
        r['actual'] = {'camara': L, 'ruta': ruta, 'ruta_centro': ruta_c, 'k': 0,
                       'n': int(comp.par.Npares), 'descartados': 0, 'nodos': 0, 'seguidos': 0,
                       'inicio': datetime.datetime.now()}
        # El seguimiento arranca desde el centro marcado y después sigue a la partícula.
        comp.par.Referenciax, comp.par.Referenciay = comp.par.Centrox.eval(), comp.par.Centroy.eval()
        comp.op('por_pasada').module.limpiar()
        self._encolar(0)

    def _escribir(self, k):
        comp, r, c = self.ownerComp, self.Corrida, self.Corrida['actual']
        campo = comp.op('bucle').op('validar').numpyArray(delayed=False)
        win = int(comp.par.Ventana)
        paso = win // 2
        ny, nx = campo.shape[:2]
        Y, X = np.meshgrid(np.arange(ny) * paso + (win - 1) * 0.5,
                           np.arange(nx) * paso + (win - 1) * 0.5, indexing='ij')
        dentro = comp.op('piv_numerico').module.dentro_roi(
            ny, nx, win, float(comp.par.Centrox), float(comp.par.Centroy), float(comp.par.Roiradio))
        s = float(comp.par.Escala) or 1.0
        dtp, tc = self._tiempo_par(k)
        m = dentro.ravel()
        x, y = X.ravel()[m], Y.ravel()[m]
        u = campo[..., 0].ravel()[m] * s / dtp
        v = campo[..., 1].ravel()[m] * s / dtp
        malo = campo[..., 2].ravel()[m] > 0.5
        c['descartados'] += int(malo.sum())
        c['nodos'] = int(m.sum())
        # Donde la validación descartó el vector no se inventa nada: va NaN. El
        # relleno con los vecinos solo sale si se pidió expresamente.
        if not comp.par.Rellenar:
            u, v = np.where(malo, np.nan, u), np.where(malo, np.nan, v)
        bloque = np.column_stack([np.full(x.size, k), np.full(x.size, tc), x, y, x * s, y * s, u, v, ~malo])
        np.savetxt(r['csv'], bloque, delimiter=',',
                   fmt=['%d', '%.5f', '%.1f', '%.1f', '%.3f', '%.3f', '%.4f', '%.4f', '%d'])

        centro = comp.op('centro')
        xc, yc, seguido = float(centro['x']), float(centro['y']), float(centro['seguido']) > 0.5
        r['csv_centro'].write('%d,%.5f,%.3f,%.3f,%.4f,%.4f,%d,%d\n' % (
            k, tc, xc, yc, xc * s, yc * s, seguido, int(centro['area'])))
        if seguido:
            c['seguidos'] += 1
            comp.par.Referenciax, comp.par.Referenciay = xc, yc     # el próximo par busca cerca

    def _siguiente(self):
        r = self.Corrida
        c = r['actual']
        c['k'] += 1
        if c['k'] >= c['n']:
            self._cerrar_camara()
            if len(r['hechas']) < len(r['camaras']):
                self._empezar_camara()
            else:
                self._cerrar()
            return
        self._encolar(c['k'])
        if c['k'] % 10 == 0:
            seg = (datetime.datetime.now() - c['inicio']).total_seconds()
            falta = seg / c['k'] * (c['n'] - c['k'])
            seguimiento = ('   partícula seguida en %d de %d' % (c['seguidos'], c['k'])
                           if self.ownerComp.par.Seguirparticula else '')
            self._say('ANALIZANDO cámara %s (%d de %d)\npar %d de %d — %.0f %%, faltan ~%.0f min%s' % (
                c['camara'], len(r['hechas']) + 1, len(r['camaras']), c['k'], c['n'],
                100.0 * c['k'] / c['n'], falta / 60, seguimiento))

    def _cerrar_camara(self):
        comp, r = self.ownerComp, self.Corrida
        c = r.get('actual')
        if c is None or r['csv'] is None:
            return
        r['csv'].close()
        r['csv_centro'].close()
        r['csv'] = r['csv_centro'] = None
        L = c['camara']
        P = comp.par
        rgb = lambda nombre: [round(float(P[L + nombre + k]), 4) for k in 'rgb']
        r['hechas'][L] = {
            'csv': os.path.basename(c['ruta']),
            'csv_centro': os.path.basename(c['ruta_centro']),
            'vista': 'superior (x, y)' if L == 'A' else 'lateral (z)',
            'pares_analizados': c['k'], 'completo': c['k'] >= c['n'],
            'calibracion': {'mm_por_px': float(P[L + 'escala']) or None,
                            'punto_a_px': [float(P[L + 'puntoax']), float(P[L + 'puntoay'])],
                            'punto_b_px': [float(P[L + 'puntobx']), float(P[L + 'puntoby'])],
                            'distancia_mm': float(P[L + 'distancia']),
                            'centro_vortice_px': [float(P[L + 'centrox']), float(P[L + 'centroy'])],
                            'radio_roi_px': float(P[L + 'roiradio']),
                            'resolucion_analisis_px': [int(comp.op('cuadro_a').width),
                                                       int(comp.op('cuadro_a').height)]},
            'colores': {'fondo': rgb('colorfondo'), 'glitter': rgb('colorglitter'),
                        'particula': rgb('colorparticula')},
            'seguimiento_particula': {'activo': bool(P.Seguirparticula), 'umbral': float(P.Umbral),
                                      'area_minima_px': int(P.Areaminima),
                                      'pares_con_particula': c['seguidos']},
            'tiempo': {'dt_por_cuadro_s': float(P.Dt), 'fps': 1.0 / float(P.Dt),
                       'cuadros_perdidos': int((self.Nint - 1).sum()),
                       'reloj_exacto': bool(self.RelojExacto),
                       'paso_temporal': ('hora de cada cuadro (reloj exacto de la simulacion)' if self.RelojExacto
                                         else 'fijo: dt por cuadro x cuadros del intervalo')},
            'vectores_por_par': c['nodos'],
            'descartados_pct': 100.0 * c['descartados'] / max(1, c['k'] * c['nodos']),
        }
        r['actual'] = None

    def _cerrar(self):
        comp, r = self.ownerComp, self.Corrida
        self.Corrida = None
        self.Cola = []
        comp.par.Corriendo = False
        self.Clave = None            # el par que se está mirando se recalcula con sus vistas
        win = int(comp.par.Ventana)
        params = {
            'grabacion': self._grabacion(),
            'fecha_analisis': r['inicio'].strftime('%Y-%m-%d %H:%M:%S'),
            'detenido_a_mano': bool(r.get('detenido')),
            'camaras': r['hechas'],
            'piv': {'metodo': 'FFT multipasada con deformacion de ventana (como PIVlab)',
                    'imagen': 'proyeccion de color sobre la direccion del glitter (sin la particula)',
                    'separacion_cuadros': self._separacion(),
                    'ventanas_px': [self._ventana(k) for k in range(self._npasadas())], 'solapamiento': 0.5,
                    'preprocesamiento': ('CLAHE' if comp.par.Imagen == 'clahe' else
                                         'lineal (proyeccion de color sin cuantizar ni recortar)'),
                    'subpixel': 'gaussiano 3 puntos', 'deformacion': 'Lanczos 8x8',
                    'validacion': 'mediana local umbral 3, desvio estandar n=8',
                    'relacion_picos_minima': float(comp.par.Picomin),
                    'rellenar_huecos': bool(comp.par.Rellenar),
                    'huecos': ('rellenados con el promedio de los vecinos' if comp.par.Rellenar
                               else 'no se rellenan: lo descartado va NaN')},
            'grabacion_info': self.Info,
            'columnas_csv': {
                'frame': 'primer cuadro del par (el segundo es frame + piv.separacion_cuadros)',
                't_s': 'centro del par, desde el primer cuadro de esa camara',
                'x_px, y_px': 'centro de la ventana; origen abajo-izquierda, y hacia arriba',
                'x_mm, y_mm': 'idem en mm (en la camara B, y es la altura z)',
                'u_mm_s, v_mm_s': 'velocidad (px/s si esa camara no tiene calibracion); '
                                  'NaN = la validacion descarto ese vector',
                'valido': '1 = medido; 0 = descartado por la validacion'},
            'columnas_csv_centro': {
                'xc_px, yc_px': 'centro del vortice en el instante medio del par (promedio de los dos cuadros)',
                'xc_mm, yc_mm': 'idem en mm',
                'seguido': '1 = particula encontrada en los dos cuadros; 0 = se repite la ultima posicion '
                           '(o el centro marcado si el seguimiento esta apagado)',
                'area_px': 'area de la mancha de la particula en el primer cuadro'},
        }
        with open(os.path.join(r['salida'], 'parametros.json'), 'w', encoding='utf-8') as f:
            json.dump(params, f, indent=2, ensure_ascii=False)
        self.Clave = None
        lineas = ['%s: %d pares%s' % (L, h['pares_analizados'], '' if h['completo'] else ' (incompleto)')
                  for L, h in r['hechas'].items()]
        self._say('%s — guardado en data glitter/%s/%s%s\n%s' % (
            'DETENIDO' if r.get('detenido') else 'LISTO', r['base'], os.path.basename(r['salida']),
            r.get('etapas', ''), '\n'.join(lineas)))

    # ------------------------------------------------------ póster
    def GuardarEtapas(self):
        """Imágenes de cada etapa del par que se ve, en el último análisis de esta grabación."""
        if self.Corrida is not None:
            self._say('Durante ANALIZAR TODO no: esperá a que termine o detenelo.')
            return
        if (self.T is None or self.Cola or self.EnCurso is not None
                or self.ownerComp.op('por_pasada').module.leer('suave', self._npasadas() - 1) is None):
            self._say('Esperá a que termine de calcularse el par y volvé a apretar.')
            return
        carpeta = self._ultima_carpeta_analisis() or self._nueva_carpeta_analisis()
        destino, n = self._exportar_etapas(carpeta)
        self._say('Imágenes guardadas (%d) en\n%s' % (
            n, os.path.relpath(destino, os.path.join(self._salida(), '..'))))

    def _exportar_etapas(self, carpeta):
        """Escribe imagenes/camara<L>_par<k>/ dentro de la carpeta del análisis. Por iteración:

          ventanas_iterN.png        zoom: cuadros A (rojo claro) y B (azul claro) sin deformar y
                                    una grilla: las ventanas de la iteración N (su tamaño real)
                                    donde se leen en el cuadro B, movidas por el predictor, y
                                    una flecha por ventana con el desplazamiento promedio medido.
          correlacion_iterN_*.png   ventana A, ventana B y plano de correlación de una ventana.
          desplazamiento_iterN.png  el desplazamiento medido en la iteración N sobre todo el ROI,
                                    con la misma escala de flechas en todas

        Solo imágenes, sin texto. Los números para rotular van en leeme.txt.
        """
        comp, P = self.ownerComp, self.ownerComp.par
        pp, dib = comp.op('por_pasada').module, comp.op('dibujo').module
        L, k = P.Camara.eval(), int(P.Cuadro)
        n = self._npasadas()
        destino = os.path.join(carpeta, 'imagenes', 'camara%s_cuadros%05d_%05d' % (L, k, k + self._separacion()))
        os.makedirs(destino, exist_ok=True)
        for viejo in os.listdir(destino):               # que no queden imágenes de otra configuración
            if viejo.endswith('.png'):
                os.remove(os.path.join(destino, viejo))
        cx, cy, radio = float(P.Centrox), float(P.Centroy), float(P.Roiradio)
        fondo = pp.leer('fondo', 0)
        H, W = fondo.shape
        A, B = pp.leer('deformada_a', 0), pp.leer('deformada_b', 0)    # iteración 1: sin deformar
        dentro_roi = comp.op('piv_numerico').module.dentro_roi

        # Región del detalle: un cuadrado de Zoompx alrededor del punto de ejemplo.
        lado = int(P.Zoompx)
        ex_, ey_ = float(P.Ejemplox), float(P.Ejemploy)
        if ex_ <= 0 and ey_ <= 0:
            ex_, ey_ = cx + 0.5 * radio, cy
        zoom = dict(x0=int(np.clip(round(ex_ - lado / 2), 0, W - lado)),
                    y0=int(np.clip(round(ey_ - lado / 2), 0, H - lado)), w=lado, h=lado,
                    escala=LADO_IMAGEN / lado)
        enc = dib.encuadre_roi(W, H, cx, cy, radio, LADO_IMAGEN)

        # Escala de flechas común: el percentil 95 del campo final mide 1.6 pasos finales (las
        # flechas de la última iteración se encadenan a lo largo de las líneas de corriente).
        final = pp.leer('campo', n - 1)
        ny, nx = final[0].shape[:2]
        medidos = dentro_roi(ny, nx, final[1], cx, cy, radio) & (final[0][..., 2] <= 0.5)
        mod = np.hypot(final[0][..., 0], final[0][..., 1])[medidos]
        p95 = float(np.percentile(mod, 95)) if mod.size else 1.0
        px_por_px = 1.6 * final[2] / max(p95, 1e-6)

        filas, notas = [], []
        corr = comp.op('dibujar_correlacion_callbacks').module
        ampliar = dib.factor_flechas(final[0], final[1])    # mismo factor en todas las iteraciones
        for pasada in range(n):
            it = pasada + 1
            campo = pp.leer('campo', pasada)
            previo = pp.leer('suave', pasada - 1) if pasada > 0 else None
            # La ventana de ejemplo de esta iteración: se guarda su correlación.
            ej = corr.ventana_de_ejemplo(comp, pp, pasada)
            destacada = None
            if ej is not None:
                wa, wb, win, du, dv, xn, yn = ej
                paso = campo[2]
                destacada = (int(round((xn - (win - 1) * 0.5) / paso)), int(round((yn - (win - 1) * 0.5) / paso)), paso)
                va, vb, plano = dib.correlacion(wa, wb, win, du, dv, LADO_VENTANA)
                dib.guardar_png(os.path.join(destino, 'correlacion_iter%d_ventanaA.png' % it), va)
                dib.guardar_png(os.path.join(destino, 'correlacion_iter%d_ventanaB.png' % it), vb)
                dib.guardar_png(os.path.join(destino, 'correlacion_iter%d_plano.png' % it), plano)
                i, j = destacada[0], destacada[1]
                notas.append('  iteracion %d: ventana de %d px centrada en (%.1f, %.1f) px; residuo medido por la '
                             'correlacion (%+.3f, %+.3f) px; desplazamiento total de la ventana (%+.3f, %+.3f) px' % (
                                 it, win, xn, yn, du, dv, campo[0][j, i, 0], campo[0][j, i, 1]))
            dib.guardar_png(os.path.join(destino, 'ventanas_iter%d.png' % it),
                            dib.ventanas(A, B, previo, self._ventana(pasada), 0, 0, 0, zoom,
                                         medido=campo, ampliar=ampliar))
            dib.guardar_png(os.path.join(destino, 'desplazamiento_iter%d.png' % (pasada + 1)),
                            dib.desplazamiento(fondo, campo[0], campo[1], campo[2], cx, cy, radio,
                                               dentro_roi, enc, px_por_px,
                                               rellenados=bool(P.Rellenar)))
            m = dentro_roi(campo[0].shape[0], campo[0].shape[1], campo[1], cx, cy, radio)
            ok = m & (campo[0][..., 2] <= 0.5)
            filas.append('  iteracion %d: ventana %3d px, paso %2d px, %4d nodos, %4d medidos, %3d descartados, '
                         'desplazamiento mediano %.2f px' % (
                             pasada + 1, campo[1], campo[2], int(m.sum()), int(ok.sum()), int(m.sum() - ok.sum()),
                             float(np.median(np.hypot(campo[0][..., 0], campo[0][..., 1])[ok])) if ok.any() else 0))

        mm = float(P.Escala)
        dtp = self._tiempo_par(k)[0]
        with open(os.path.join(destino, 'leeme.txt'), 'w', encoding='utf-8') as f:
            f.write('Imagenes del PIV de un par de cuadros (sin rotulos)\n')
            f.write('grabacion %s, camara %s, par %d (cuadros %d y %d), dt del par %.5f s\n\n' % (
                os.path.basename(self._grabacion()), L, k, k, k + self._separacion(), dtp))
            f.write('ventanas_iterN.png  (lo que considera la iteracion N)\n'
                    '  Region de %d x %d px de la imagen (x = %d..%d, y = %d..%d, y hacia arriba), ampliada a %d px (x%.2f).\n'
                    '  Fondo: cuadro A en rojo claro y cuadro B en azul claro, SIN deformar (gris = coinciden):\n'
                    '  cada destello aparece dos veces, separado por lo que se movio entre los dos cuadros.\n'
                    '  Grilla negra: las ventanas de interrogacion de la iteracion, de su tamano real (una si y una\n'
                    '  no: con 50 %% de solapamiento hay otra ventana centrada en cada esquina), en el lugar donde se\n'
                    '  leen en el cuadro B: p + d(p)/2, con d = predictor de la iteracion = campo (suavizado) de la\n'
                    '  iteracion anterior, a escala real. En el cuadro A se leen en p - d(p)/2 (la misma grilla\n'
                    '  corrida hacia el otro lado); la correlacion compara las dos lecturas.\n'
                    '  iter1: sin predictor, la grilla es recta. En las siguientes la grilla se mueve con el flujo;\n'
                    '  debajo, en gris claro, la misma grilla sin mover como referencia.\n'
                    '  Flechas violetas: el desplazamiento promedio que calculo la iteracion N en cada ventana\n'
                    '  dibujada, centradas en la ventana y AMPLIADAS x%g (igual en todas las iteraciones): apuntan\n'
                    '  hacia donde se movieron los destellos de esa ventana, de A (rojo) a B (azul).\n'
                    '' % (
                        lado, lado, zoom['x0'], zoom['x0'] + lado, zoom['y0'], zoom['y0'] + lado, LADO_IMAGEN,
                        zoom['escala'], ampliar))
            if mm > 0:
                f.write('  Escala de ventanas_*: 1 px = %.4f mm (%.1f px por cm).\n' % (
                    mm / zoom['escala'], 10.0 * zoom['escala'] / mm))
            f.write('\n')
            f.write('correlacion_iterN_ventanaA.png, correlacion_iterN_ventanaB.png, correlacion_iterN_plano.png\n'
                    '  La ventana de ejemplo (posicion abajo), tal como la compara la iteracion N: leida en el cuadro A\n'
                    '  (ventanaA) y en el B (ventanaB), ya corridas por el predictor. Pixeles sin suavizar,\n'
                    '  %d x %d px cada imagen.\n'
                    '  plano = correlacion cruzada de las dos ventanas (calculada con FFT): cuanto se parecen\n'
                    '  B y A si se corre B cada desplazamiento posible. Centro (cruz gris) = desplazamiento 0;\n'
                    '  hacia la derecha +x, hacia arriba +y. 1 px de desplazamiento = (%d / ventana) px del plano:\n'
                    '  %s. Cuadrado blanco = zona de busqueda (|d| <= ventana/4).\n'
                    '  Anillo = pico con ajuste gaussiano de 3 puntos: el desplazamiento que falta sumarle al\n'
                    '  predictor (residuo). En la iter1 (sin predictor) el pico esta lejos del centro; en las\n'
                    '  siguientes queda cerca: el predictor ya explico casi todo el movimiento.\n' % (
                        LADO_VENTANA, LADO_VENTANA, LADO_VENTANA,
                        ', '.join('iter%d %g px' % (p + 1, LADO_VENTANA / self._ventana(p)) for p in range(n))))
            f.write('\n'.join(notas) + '\n\n')
            f.write('desplazamiento_iterN.png\n'
                    '  Todo el ROI: recorte x = %d..%d, y = %d..%d px, ampliado x%.3f (%d px).\n'
                    '  Flechas azules = desplazamiento medido en esa iteracion. Donde la validacion descarto\n'
                    '  el vector no hay flecha: ahi no se midio nada (con "Rellenar huecos" prendido aparecen\n'
                    '  en rojo, que es el promedio de los vecinos). Misma escala en todas: 1 px = %.2f px de la\n'
                    '  imagen original = %.1f px del PNG.\n' % (
                        enc['x0'], enc['x0'] + enc['w'], enc['y0'], enc['y0'] + enc['h'], enc['escala'], LADO_IMAGEN,
                        px_por_px, px_por_px * enc['escala']))
            if mm > 0:
                f.write('  Escala del PNG: 1 px = %.4f mm (%.2f px por cm). 1 px de desplazamiento en este par = %.2f mm/s.\n' % (
                    mm / enc['escala'], 10.0 * enc['escala'] / mm, mm / dtp))
            f.write('\nIteraciones:\n' + '\n'.join(filas) + '\n')
        return destino, 5 * n
