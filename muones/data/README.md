# Datos de decaimientos de muones (DRS4)

Trazas de los tres centelladores guardadas por `drs_exam` (ver `muones/grupo_7_jueves`)
entre el 11/9 y el 29/10 de 2026. Son los eventos `_selected` de cada corrida: los que
pasaron la selección en línea del programa de adquisición.

Las trazas **no están en git** (31 GB): están en el release
[`muones-decaimientos-v1`](https://github.com/olistarosta/labo-5-2026/releases/tag/muones-decaimientos-v1).
En git solo quedan este archivo, `MANIFEST.csv` y dos archivos sueltos de ejemplo
(`run_20260911_135116_idx0_selected.txt`, `run_20261022_194538_idx0_selected.txt`).

## Bajar los datos

Desde la raíz del repositorio:

```bash
gh release download muones-decaimientos-v1 -R olistarosta/labo-5-2026 -D muones/data -p "run_*.zip"
```

y después descomprimir todos los zips **dentro de esta carpeta** (cada zip ya trae su
carpeta `run_.../`). En PowerShell:

```powershell
Get-ChildItem muones/data/run_*.zip | ForEach-Object { Expand-Archive $_ muones/data }
```

El release también trae `indice.zip`: el índice de pulsos ya calculado (61 MB), que se
descomprime acá mismo como `_indice/`. Con eso el notebook
`muones/analisis/decaimientos.ipynb` no necesita recorrer los 31 GB la primera vez.

## Estructura

```
data/
  MANIFEST.csv                          un renglón por archivo (bytes, CRC32, carpeta de origen en Drive)
  run_<AAAAMMDD>_<HHMMSS>/              una carpeta por corrida (una ejecución de drs_exam)
    run_<...>_idx<N>_selected.txt       un archivo cada ~2000 disparos
  _indice/<run>.npz                     caché que arma decaimientos.py (no va a git)
```

| corrida | archivos | eventos | primer evento | último evento | tamaño |
|---|---:|---:|---|---|---:|
| `run_20260911_135116` | 1728 | 235 402 | 11/9 13:51 | 17/9 14:11 | 7,3 GB |
| `run_20260918_135753` | 1425 | 103 809 | 18/9 13:57 | 21/9 22:25 | 3,3 GB |
| `run_20260923_175832` | 241 | 34 597 | 23/9 17:58 | 24/9 14:03 | 1,1 GB |
| `run_20260925_135327` | 1635 | 196 586 | 25/9 21:09 | 1/10 14:07 | 6,1 GB |
| `run_20261001_141452` | 1 | 13 | 1/10 14:14 | 1/10 14:15 | 0,4 MB |
| `run_20261002_135657` | 1736 | 201 978 | 3/10 06:44 | 9/10 08:06 | 6,3 GB |
| `run_20261016_134329` | 1724 | 100 152 | 16/10 13:43 | 22/10 14:03 | 3,2 GB |
| `run_20261023_134416` | 1 | 7 | 23/10 13:44 | 23/10 13:45 | 0,2 MB |
| `run_20261023_134533` | 1725 | 105 796 | 23/10 13:45 | 29/10 14:07 | 3,3 GB |
| **total** | **10 216** | **978 340** | | | **31 GB** |

Cosas a saber:

- `run_20260918_135753`: los archivos `idx963` a `idx1424` (462) están vacíos (0 bytes);
  la corrida dejó de escribir el 21/9 a las 22:25.
- `run_20260925_135327` empieza en `idx87` y `run_20261002_135657` en `idx201`: los
  primeros archivos de esas corridas no estaban en la carpeta de Drive. Tampoco estaban
  `idx88`, `89`, `92`, `94` y `1503` de la primera, ni `idx202`–`204` y `1369` de la segunda.
- Las dos corridas de un solo archivo (`run_20261001_141452`, `run_20261023_134416`) son
  pruebas de menos de un minuto.

## Formato de los archivos

Texto, un bloque por evento:

```
# DRS4 export. File index 0, run_id=./medidas/run_20260918_135753
# Columns: t[ns] u1[mV] u2[mV] u3[mV]
Event #5  ts=2026-09-18 13:57:58.446  [SELECTED]
  t[ns]  u1[mV]  u2[mV]  u3[mV]
  0.000   -17.2   -10.0   -11.7
  1.660   -15.0    -7.0   -11.2
  ...                                  (1024 muestras, ~1,7 µs a 0,6 GS/s)
```

`Event #` es el número de disparo dentro de la corrida (cuenta también los no
seleccionados) y `ts` la hora de la Raspberry. El disparo es por `u1` y el pulso que
dispara cae siempre cerca de los 207 ns.

## De dónde salieron

Bajados de la carpeta `Datos` de Drive el 9/11/2026 (12 zips). Se reorganizaron por
corrida; la columna `carpeta_drive` de `MANIFEST.csv` dice de qué subcarpeta venía cada
archivo. Se descartaron 75 duplicados exactos (mismo CRC32): 74 archivos de
`datos25sep_1oct` que estaban repetidos en `Últimos_datos`, y
`run_20261016_134329_idx1631_selected(1).txt`.
