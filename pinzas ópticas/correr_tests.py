"""Corre todos los tests con el Python de TouchDesigner:

    "C:/Program Files/Derivative/TouchDesigner/bin/python.exe" correr_tests.py [-k algo]

(pytest está en TouchDesigner/lib/, que se agrega al FINAL de sys.path para no tapar el numpy de TD).
"""
import os
import sys

RAIZ = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, RAIZ)
sys.path.append(os.path.join(RAIZ, 'TouchDesigner', 'lib'))
import pytest  # noqa: E402

sys.exit(pytest.main([os.path.join(RAIZ, 'tests'), '-q', '-p', 'no:cacheprovider'] + sys.argv[1:]))
