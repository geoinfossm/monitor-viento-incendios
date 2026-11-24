"""Pruebas de la lógica geométrica y de viento. Ejecutar:  python test_logica.py"""
import math
import monitor_viento_incendios as m

# Distancia: 1 grado de latitud ≈ 111 km
assert abs(m.haversine_km(-35, -71, -34, -71) - 111.2) < 0.5

# Rumbos
assert abs(m.rumbo_deg(-35, -71, -34, -71) - 0) < 0.5      # norte
assert abs(m.rumbo_deg(-35, -71, -35, -70) - 90) < 1       # este
assert abs(m.rumbo_deg(-35, -71, -36, -71) - 180) < 0.5    # sur

# Diferencia angular
assert m.dif_ang(350, 10) == 20
assert m.dif_ang(0, 180) == 180
assert m.dif_ang(90, 270) == 180

# Promedio de direcciones: 350° y 10° deben dar ~0°, NO 180°
r = m.combinar_modelos([(0.5, 10, 350, 15), (0.5, 10, 10, 15)])
assert min(r["dir_desde"], 360 - r["dir_desde"]) < 1, r
assert abs(r["vel_kmh"] - 10) < 1e-9

# Modelos que se anulan -> dirección indefinida
r = m.combinar_modelos([(0.5, 20, 90, None), (0.5, 20, 270, None)])
assert r["dir_desde"] is None

# Un modelo sin dato redistribuye el peso
r = m.combinar_modelos([(0.4, 10, 180, 20), (0.3, None, None, None), (0.3, 20, 180, 30)])
assert abs(r["vel_kmh"] - (0.4 * 10 + 0.3 * 20) / 0.7) < 1e-9

crit = m.CFG["criterios"]

# Fuego al NORTE del establecimiento (rumbo 0). Viento DESDE el norte (0°) sopla hacia el sur,
# es decir, hacia el establecimiento -> delta 0 -> alineado
assert m.clasificar(8, m.dif_ang(0, 0), 20, crit) == "ALERTA"
assert m.clasificar(20, m.dif_ang(0, 0), 20, crit) == "VIGILANCIA"
assert m.clasificar(28, m.dif_ang(0, 0), 20, crit) == "INFORMATIVO"
# Viento DESDE el sur (180°): empuja el fuego lejos del establecimiento
assert m.clasificar(8, m.dif_ang(180, 0), 20, crit) == "INFORMATIVO"
# Muy cerca, aunque el viento lo aleje
assert m.clasificar(3, m.dif_ang(180, 0), 20, crit) == "VIGILANCIA"
# Viento débil no escala por alineación
assert m.clasificar(8, 0, 3, crit) == "INFORMATIVO"

assert m.tipo_viento(10, 45) == "HACIA_ESTABLECIMIENTO"
assert m.tipo_viento(90, 45) == "LATERAL"
assert m.tipo_viento(170, 45) == "ALEJANDOSE"
assert m.cardinal(0) == "N" and m.cardinal(270) == "O" and m.cardinal(225) == "SO"

print("OK: todas las pruebas de lógica pasaron")
