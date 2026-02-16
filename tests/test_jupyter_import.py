#!/usr/bin/env python3
"""Test pour diagnostiquer le problème d'import dans Jupyter"""

print("=== Test d'importation EPyR ===")

# Test 1: Import standard
try:
    import epyr

    print(f"✅ Import epyr réussi - Version: {epyr.__version__}")
except Exception as e:
    print(f"❌ Erreur import epyr: {e}")

# Test 2: Vérification des fonctions
functions_to_check = ["plot_1d", "plot_2d_map", "plot_2d_waterfall"]
for func_name in functions_to_check:
    available = hasattr(epyr, func_name)
    print(f"   {func_name}: {'✅' if available else '❌'}")

# Test 3: Import direct du module eprplot
try:
    from epyr import eprplot

    print("✅ Import direct de epyr.eprplot réussi")
    print(f"   Fonctions dans eprplot: {eprplot.__all__}")
except Exception as e:
    print(f"❌ Erreur import epyr.eprplot: {e}")

# Test 4: Test des fonctions
try:
    import numpy as np

    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    fig, ax = epyr.plot_1d(x, y, {"XAXIS_NAME": "Test", "XAXIS_UNIT": "a.u."})
    print("✅ Fonction plot_1d fonctionne")

    # Fermer la figure pour éviter l'affichage
    import matplotlib.pyplot as plt

    plt.close(fig)

except Exception as e:
    print(f"❌ Erreur test plot_1d: {e}")

print("\n=== Instructions pour Jupyter ===")
print("Si ce script fonctionne mais pas Jupyter, essayez dans Jupyter:")
print("1. Kernel -> Restart & Clear Output")
print("2. Ou ajoutez en première cellule:")
print("   import sys")
print("   print(sys.executable)")
print("   print(sys.path[:3])")
print("3. Ou essayez:")
print("   import importlib")
print("   import epyr")
print("   importlib.reload(epyr)")
