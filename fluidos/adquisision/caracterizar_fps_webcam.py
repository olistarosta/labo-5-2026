"""
Caracterización de FPS de una webcam.

Captura frames de una webcam, registra el timestamp (tiempo de la
computadora) en el que se recibe cada frame, calcula los intervalos
entre frames consecutivos, los grafica en un histograma y guarda
los datos en un CSV.

Uso:
    python caracterizar_fps_webcam.py --duracion 10 --camara 0
"""

import argparse
import csv
import time

import cv2
import matplotlib.pyplot as plt
import numpy as np


def capturar_timestamps(camara_idx: int, duracion_s: float) -> list[float]:
    """Captura frames de la webcam durante 'duracion_s' segundos y
    devuelve la lista de timestamps (time.perf_counter()) en los que
    se recibió cada frame."""

    cap = cv2.VideoCapture(camara_idx)
    if not cap.isOpened():
        raise RuntimeError(f"No se pudo abrir la cámara {camara_idx}")

    timestamps = []
    t_inicio = time.perf_counter()

    try:
        while (time.perf_counter() - t_inicio) < duracion_s:
            ret, frame = cap.read()
            t = time.perf_counter()  # timestamp apenas llega el frame
            if not ret:
                print("Advertencia: frame no recibido, se omite.")
                continue
            timestamps.append(t)
    finally:
        cap.release()

    return timestamps


def calcular_intervalos(timestamps: list[float]) -> np.ndarray:
    """Calcula los intervalos (deltas) entre timestamps consecutivos, en ms."""
    ts = np.array(timestamps)
    intervalos_s = np.diff(ts)
    return intervalos_s * 1000.0  # a milisegundos


def guardar_csv(timestamps: list[float], intervalos_ms: np.ndarray, path_csv: str) -> None:
    """Guarda timestamps e intervalos en un CSV.

    Columnas: frame_idx, timestamp_s, intervalo_ms (vacío para el primer frame)
    """
    with open(path_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["frame_idx", "timestamp_s", "intervalo_ms"])
        for i, t in enumerate(timestamps):
            intervalo = intervalos_ms[i - 1] if i > 0 else ""
            writer.writerow([i, f"{t:.6f}", intervalo])

    print(f"CSV guardado en: {path_csv}")


def graficar_histograma(intervalos_ms: np.ndarray, path_png: str) -> None:
    """Grafica un histograma de los intervalos entre frames y lo guarda como PNG."""
    fps_promedio = 1000.0 / np.mean(intervalos_ms)

    plt.figure(figsize=(8, 5))
    plt.hist(intervalos_ms, bins=40, color="steelblue", edgecolor="black")
    plt.axvline(np.mean(intervalos_ms), color="red", linestyle="--",
                label=f"Media: {np.mean(intervalos_ms):.2f} ms")
    plt.xlabel("Intervalo entre frames (ms)")
    plt.ylabel("Frecuencia")
    plt.title(f"Histograma de intervalos entre frames\nFPS promedio ≈ {fps_promedio:.2f}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path_png, dpi=150)
    print(f"Histograma guardado en: {path_png}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Caracteriza el FPS de una webcam")
    parser.add_argument("--camara", type=int, default=0, help="Índice de la cámara (default: 0)")
    parser.add_argument("--duracion", type=float, default=10.0, help="Duración de captura en segundos")
    parser.add_argument("--csv", type=str, default="fps_webcam.csv", help="Path del CSV de salida")
    parser.add_argument("--png", type=str, default="histograma_fps.png", help="Path del histograma de salida")
    args = parser.parse_args()

    print(f"Capturando de la cámara {args.camara} durante {args.duracion} s...")
    timestamps = capturar_timestamps(args.camara, args.duracion)

    n_frames = len(timestamps)
    if n_frames < 2:
        raise RuntimeError("No se capturaron suficientes frames para calcular intervalos.")

    intervalos_ms = calcular_intervalos(timestamps)

    print(f"Frames capturados: {n_frames}")
    print(f"Intervalo medio: {np.mean(intervalos_ms):.3f} ms  "
          f"(FPS promedio ≈ {1000.0/np.mean(intervalos_ms):.2f})")
    print(f"Desvío estándar: {np.std(intervalos_ms):.3f} ms")
    print(f"Mínimo: {np.min(intervalos_ms):.3f} ms | Máximo: {np.max(intervalos_ms):.3f} ms")

    guardar_csv(timestamps, intervalos_ms, args.csv)
    graficar_histograma(intervalos_ms, args.png)


if __name__ == "__main__":
    main()
