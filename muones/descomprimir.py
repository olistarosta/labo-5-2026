#!/usr/bin/env python3
"""
Descomprime los archivos comprimidos que haya en la carpeta data (junto a este script).
Soporta .zip, .tar, .tar.gz, .tgz, .tar.bz2. Cada archivo se extrae en una carpeta con su mismo nombre.

Uso (desde la carpeta muones):
    python3 descomprimir.py
    python3 descomprimir.py /otra/carpeta
"""
import os
import sys
import tarfile
import zipfile

CARPETA_SCRIPT = os.path.dirname(os.path.abspath(__file__))
carpeta = sys.argv[1] if len(sys.argv) > 1 else os.path.join(CARPETA_SCRIPT, "data")

print(f"Carpeta buscada: {carpeta}")
if not os.path.isdir(carpeta):
    sys.exit("ERROR: esa carpeta no existe. Creala (muones/data) y poné ahí los archivos comprimidos.")

TAR_EXT = (".tar.gz", ".tgz", ".tar.bz2", ".tar")
ZIP_EXT = (".zip",)
NO_SOPORTADOS = (".7z", ".rar")

todos = []
for raiz, _, nombres in os.walk(carpeta):
    for n in nombres:
        todos.append(os.path.join(raiz, n))

print(f"Archivos encontrados en total: {len(todos)}")
if not todos:
    sys.exit("La carpeta está vacía: no hay nada para descomprimir.")

hechos = 0
for ruta in sorted(todos):
    bajo = ruta.lower()
    destino = None
    try:
        if bajo.endswith(ZIP_EXT):
            destino = ruta[:-4]
            with zipfile.ZipFile(ruta) as z:
                z.extractall(destino)
        elif bajo.endswith(TAR_EXT):
            for ext in TAR_EXT:
                if bajo.endswith(ext):
                    destino = ruta[:-len(ext)]
                    break
            with tarfile.open(ruta) as t:
                t.extractall(destino)
        elif bajo.endswith(NO_SOPORTADOS):
            print(f"  NO SOPORTADO (usá The Unarchiver o 7-Zip): {ruta}")
            continue
        else:
            continue
    except Exception as e:
        print(f"  ERROR al descomprimir {ruta}: {e}")
        continue
    hechos += 1
    print(f"  descomprimido: {ruta}  ->  {destino}")

if hechos == 0:
    print("\nNo se descomprimió nada: ningún archivo tiene extensión .zip, .tar, .tar.gz, .tgz o .tar.bz2.")
    print("Estos son los primeros archivos que vi (fijate la extensión):")
    for r in sorted(todos)[:10]:
        print("   ", os.path.relpath(r, carpeta))
else:
    print(f"\nListo: {hechos} archivo(s) descomprimido(s).")