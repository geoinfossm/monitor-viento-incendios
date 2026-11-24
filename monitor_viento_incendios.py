#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Monitor viento–incendios por establecimiento de salud (UGREyD, SS Maule)

Qué hace en cada corrida:
  1. Lee la tabla de establecimientos (data/establecimientos.csv).
  2. Lee los focos de incendio de CONAF (URL configurable).
  3. Obtiene viento (velocidad, ráfaga, dirección) con el ensamble ECMWF/GFS/ICON
     vía Open-Meteo, en celdas de ~10 km (para no repetir consultas).
  4. Para cada establecimiento calcula, por cada incendio dentro del radio de
     búsqueda: distancia, rumbo hacia el fuego y si el viento sopla desde el
     fuego hacia el establecimiento.
  5. Publica capas GeoJSON listas para ArcGIS Online (carpeta site/).

Herramienta de tamizaje: los modelos tienen resolución de 9–25 km y el viento
local (valles, cordillera, pendiente) puede diferir. Sirve para priorizar a qué
establecimientos mirar primero, no para predecir la trayectoria del fuego.
"""
import csv
import json
import math
import os
import sys
import time
import unicodedata
from datetime import datetime, timezone, timedelta
from pathlib import Path

import requests

try:
    from zoneinfo import ZoneInfo
    TZ_CL = ZoneInfo("America/Santiago")
except Exception:  # pragma: no cover
    TZ_CL = timezone(timedelta(hours=-3))

AQUI = Path(__file__).resolve().parent
CFG = json.loads((AQUI / "config.json").read_text(encoding="utf-8"))
MOCK = os.environ.get("MOCK") == "1"
SALIDA = Path(os.environ.get("SALIDA_DIR", AQUI / "site"))
CACHE_DIR = Path(os.environ.get("CACHE_DIR", AQUI / "cache"))
OPEN_METEO = "https://api.open-meteo.com/v1/forecast"
ADVERTENCIAS = []

RANK = {"ALERTA": 3, "VIGILANCIA": 2, "INFORMATIVO": 1, "SIN_INCENDIO": 0}
CARDINALES = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
              "S", "SSO", "SO", "OSO", "O", "ONO", "NO", "NNO"]


def log(msg):
    print(msg, flush=True)


def advertir(msg):
    ADVERTENCIAS.append(msg)
    log("[!] " + msg)


# ---------------------------------------------------------------------------
# Geometría y viento (funciones puras, probadas en test_logica.py)
# ---------------------------------------------------------------------------
def haversine_km(la1, lo1, la2, lo2):
    r = 6371.0088
    p1, p2 = math.radians(la1), math.radians(la2)
    dphi = p2 - p1
    dl = math.radians(lo2 - lo1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def rumbo_deg(la1, lo1, la2, lo2):
    """Rumbo inicial (0=N, 90=E) desde el punto 1 hacia el punto 2."""
    p1, p2 = math.radians(la1), math.radians(la2)
    dl = math.radians(lo2 - lo1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return (math.degrees(math.atan2(y, x)) + 360) % 360


def dif_ang(a, b):
    """Diferencia angular mínima entre dos rumbos, 0–180."""
    d = abs(a - b) % 360
    return d if d <= 180 else 360 - d


def cardinal(deg):
    if deg is None:
        return ""
    return CARDINALES[int((deg + 11.25) // 22.5) % 16]


def combinar_modelos(obs):
    """
    obs: lista de (peso, vel_kmh, dir_desde_deg, rafaga_kmh) por modelo.
    La dirección NO se promedia como número (350° y 10° darían 180°): se
    ponderan las componentes u/v y se reconstruye el ángulo.
    La velocidad mostrada es el promedio ponderado de las velocidades.
    """
    obs = [o for o in obs if o[1] is not None and o[2] is not None]
    if not obs:
        return None
    den = sum(o[0] for o in obs)
    if den <= 0:
        return None
    u = sum(o[0] * (-o[1] * math.sin(math.radians(o[2]))) for o in obs) / den
    v = sum(o[0] * (-o[1] * math.cos(math.radians(o[2]))) for o in obs) / den
    vel = sum(o[0] * o[1] for o in obs) / den
    dir_desde = None if math.hypot(u, v) < 0.5 else (math.degrees(math.atan2(-u, -v)) + 360) % 360
    rafs = [(o[0], o[3]) for o in obs if o[3] is not None]
    raf = (sum(w * g for w, g in rafs) / sum(w for w, _ in rafs)) if rafs else None
    disp_dir = None
    if dir_desde is not None and len(obs) >= 2:
        disp_dir = max(dif_ang(o[2], dir_desde) for o in obs)
    disp_vel = (max(o[1] for o in obs) - min(o[1] for o in obs)) if len(obs) >= 2 else None
    return {
        "vel_kmh": vel, "rafaga_kmh": raf, "dir_desde": dir_desde,
        "dir_hacia": None if dir_desde is None else (dir_desde + 180) % 360,
        "disp_dir": disp_dir, "disp_vel": disp_vel, "n_modelos": len(obs),
    }


def tipo_viento(delta, tol):
    if delta is None:
        return "INDEFINIDO"
    if delta <= tol:
        return "HACIA_ESTABLECIMIENTO"
    if delta >= 180 - tol:
        return "ALEJANDOSE"
    return "LATERAL"


def clasificar(dist, delta, vel, crit):
    """Clasificación por incendio. Umbrales en config.json (a validar con COGRID/CONAF)."""
    alineado = delta is not None and delta <= crit["tolerancia_grados"]
    if alineado and vel is not None and vel >= crit["viento_minimo_kmh"]:
        if dist <= crit["radio_alerta_km"]:
            return "ALERTA"
        if dist <= crit["radio_vigilancia_km"]:
            return "VIGILANCIA"
    if dist <= crit["radio_proximidad_km"]:
        return "VIGILANCIA"
    return "INFORMATIVO"


# ---------------------------------------------------------------------------
# Utilidades de lectura
# ---------------------------------------------------------------------------
def _norm(s):
    return unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().lower().strip()


def _num(x):
    try:
        return float(str(x).strip().replace(",", "."))
    except (ValueError, TypeError):
        return None


ALIAS = {
    "id": ["id", "codigo", "cod", "codigo_vigente", "codigo_establecimiento", "id_establecimiento", "objectid"],
    "nombre": ["nombre", "establecimiento", "nombre_establecimiento", "nombre_oficial", "name"],
    "tipo": ["tipo", "tipo_establecimiento", "nivel", "tipologia", "tipo_est"],
    "comuna": ["comuna", "nombre_comuna", "comuna_nombre"],
    "lat": ["lat", "latitud", "latitude", "y", "coord_y", "lat_decimal"],
    "lon": ["lon", "lng", "long", "longitud", "longitude", "x", "coord_x", "lon_decimal"],
}


def cargar_establecimientos():
    ruta = AQUI / CFG["establecimientos_csv"]
    if not ruta.exists():
        sys.exit(f"No se encontró {ruta}. Sube tu tabla de establecimientos con ese nombre.")
    lineas = ruta.read_text(encoding="utf-8-sig").splitlines()
    delim = ";" if lineas[0].count(";") > lineas[0].count(",") else ","
    lector = csv.DictReader(lineas, delimiter=delim)
    cab = {_norm(c): c for c in (lector.fieldnames or [])}

    def col(clave, obligatoria=False):
        pedido = _norm(CFG["columnas"].get(clave, ""))
        if pedido and pedido in cab:
            return cab[pedido]
        for a in ALIAS[clave]:
            if a in cab:
                return cab[a]
        if obligatoria:
            sys.exit(f"No encuentro la columna '{clave}'. Columnas del CSV: {list(cab.values())}. "
                     f"Ajusta 'columnas' en config.json.")
        return None

    c_id, c_nom, c_tipo, c_com = col("id"), col("nombre"), col("tipo"), col("comuna")
    c_lat, c_lon = col("lat", True), col("lon", True)

    est, corregidos, descartados = [], 0, 0
    for n, fila in enumerate(lector, start=1):
        lat, lon = _num(fila.get(c_lat)), _num(fila.get(c_lon))
        if lat is None or lon is None:
            descartados += 1
            continue
        if 17 < lat < 56:      # signo omitido
            lat, corregidos = -lat, corregidos + 1
        if 66 < lon < 76:
            lon, corregidos = -lon, corregidos + 1
        if -76 < lat < -66 and -56 < lon < -17:   # lat/lon invertidos
            lat, lon, corregidos = lon, lat, corregidos + 1
        if not (-56 < lat < -17 and -76 < lon < -66):
            descartados += 1
            continue
        est.append({
            "id": (fila.get(c_id) if c_id else None) or f"E{n}",
            "nombre": (fila.get(c_nom) if c_nom else "") or "",
            "tipo": (fila.get(c_tipo) if c_tipo else "") or "",
            "comuna": (fila.get(c_com) if c_com else "") or "",
            "lat": lat, "lon": lon,
        })
    if corregidos:
        advertir(f"{corregidos} coordenadas con signo omitido o lat/lon invertidos fueron corregidas automáticamente; revisa el CSV.")
    if descartados:
        advertir(f"{descartados} filas del CSV se descartaron por coordenadas vacías o fuera de Chile.")
    if not est:
        sys.exit("No quedaron establecimientos válidos.")
    log(f"Establecimientos cargados: {len(est)}")
    return est


# ---------------------------------------------------------------------------
# Incendios CONAF
# ---------------------------------------------------------------------------
CAND_CAMPOS = {
    "nombre": ["nombre", "name", "nom_incend", "nombre_incendio", "incendio"],
    "estado": ["estado", "status", "estado_fuego", "situacion"],
    "superficie": ["superficie", "sup_ha", "hectareas", "superficie_ha", "area_ha", "sup_total"],
    "fecha": ["fecha", "fecha_inicio", "fecha_hora", "inicio", "datetime", "f_inicio"],
}


def _prop(props, clave):
    cfg = CFG["conaf"]["campos"].get(clave, "")
    bajo = {_norm(k): k for k in props}
    if cfg and _norm(cfg) in bajo:
        return props[bajo[_norm(cfg)]]
    for c in CAND_CAMPOS[clave]:
        if c in bajo:
            return props[bajo[c]]
    return None


def _centroide(geom):
    if geom["type"] == "Point":
        return geom["coordinates"][1], geom["coordinates"][0]
    pts = []

    def rec(x):
        if isinstance(x[0], (int, float)):
            pts.append(x)
        else:
            for i in x:
                rec(i)
    rec(geom["coordinates"])
    return sum(p[1] for p in pts) / len(pts), sum(p[0] for p in pts) / len(pts)


def _fecha_iso(v):
    if isinstance(v, (int, float)) and v > 1e11:
        return datetime.fromtimestamp(v / 1000, timezone.utc).isoformat(timespec="minutes")
    return None if v is None else str(v)


def _leer_csv_incendios(url):
    """Lee un CSV de incendios (columnas lat/lon) y lo entrega como lista de 'features' tipo GeoJSON."""
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    texto = r.content.decode("utf-8-sig", errors="replace")
    lector = csv.DictReader(texto.splitlines())
    feats, descartadas = [], 0
    for fila in lector:
        fb = {_norm(k): v for k, v in fila.items() if k}
        lat = _num(fb.get("lat") or fb.get("latitud"))
        lon = _num(fb.get("lon") or fb.get("longitud") or fb.get("lng"))
        if lat is None or lon is None or not (-90 <= lat <= 90 and -180 <= lon <= 180):
            descartadas += 1
            continue
        if lon > 0:      # Chile está en longitudes negativas
            lon = -lon
        if lat > 0:
            lat = -lat
        props = dict(fila)
        if "sup_total" in fb:
            props["sup_total"] = _num(fb["sup_total"])
        feats.append({"geometry": {"type": "Point", "coordinates": [lon, lat]}, "properties": props})
    if descartadas:
        advertir(f"CONAF CSV: {descartadas} filas sin coordenadas válidas fueron descartadas.")
    if not feats:
        raise RuntimeError("El CSV no trajo filas con coordenadas (¿cambió el formato?).")
    return feats


def cargar_incendios(est):
    """Devuelve (lista_de_incendios_activos, fuente_ok)."""
    if MOCK:
        return mock_incendios(est), True
    url = os.environ.get("CONAF_URL") or CFG["conaf"]["url"]
    if not url:
        advertir("No hay URL de CONAF configurada (variable CONAF_URL o config.json).")
        return [], False
    params = None
    if ("FeatureServer" in url or "MapServer" in url) and "/query" not in url:
        url = url.rstrip("/") + "/query"
        params = {"where": CFG["conaf"].get("where", "1=1"), "outFields": "*", "outSR": 4326,
                  "f": "geojson", "returnGeometry": "true"}
    try:
        if url.lower().split("?")[0].endswith(".csv"):
            feats = _leer_csv_incendios(url)
        else:
            j = _get_json(url, params)
            if isinstance(j, dict) and "error" in j:
                raise RuntimeError(f"El servicio respondió error: {j['error']}")
            feats = j.get("features", [])
    except Exception as e:  # noqa: BLE001
        advertir(f"No se pudo leer CONAF: {e}")
        return [], False

    activos_cfg = {_norm(x) for x in CFG["conaf"]["estados_activos"]}
    estados_vistos, fuegos = {}, []
    la_min = min(e["lat"] for e in est) - 0.7
    la_max = max(e["lat"] for e in est) + 0.7
    lo_min = min(e["lon"] for e in est) - 0.7
    lo_max = max(e["lon"] for e in est) + 0.7
    for f in feats:
        g = f.get("geometry")
        if not g:
            continue
        p = f.get("properties") or {}
        lat, lon = _centroide(g)
        if not (la_min <= lat <= la_max and lo_min <= lon <= lo_max):
            continue
        estado = _prop(p, "estado")
        estados_vistos[str(estado)] = estados_vistos.get(str(estado), 0) + 1
        if estado is not None and _norm(estado) not in activos_cfg:
            continue
        fuegos.append({"nombre": _prop(p, "nombre") or "", "estado": str(estado or ""),
                       "superficie_ha": _prop(p, "superficie"), "fecha": _fecha_iso(_prop(p, "fecha")),
                       "lat": lat, "lon": lon})
    log(f"CONAF: {len(feats)} registros; en el área de estudio {sum(estados_vistos.values())}; "
        f"activos según config: {len(fuegos)}. Estados vistos: {estados_vistos}")
    if estados_vistos and not fuegos and all(_norm(k) not in activos_cfg for k in estados_vistos):
        advertir(f"Ningún estado coincide con 'estados_activos' de config.json. Estados vistos en el área: {list(estados_vistos)}")
    return fuegos, True


# ---------------------------------------------------------------------------
# Viento (Open-Meteo, ensamble por componentes u/v)
# ---------------------------------------------------------------------------
def _get_json(url, params=None, intentos=3):
    ultimo = None
    for i in range(intentos):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code in (429, 500, 502, 503, 504):
                raise RuntimeError(f"HTTP {r.status_code}")
            r.raise_for_status()
            return r.json()
        except Exception as e:  # noqa: BLE001
            ultimo = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"{ultimo}")


def clave_celda(lat, lon):
    return f"{lat:.2f},{lon:.2f}"


def celda_de(lat, lon):
    g = CFG["viento"]["celda_grados"]
    return round(round(lat / g) * g, 2), round(round(lon / g) * g, 2)


def _hora_actual(item):
    h = item.get("hourly") or {}
    tiempos = h.get("time") or []
    if not tiempos:
        return None
    ahora = datetime.now(timezone.utc)
    cand = []
    for i, t in enumerate(tiempos):
        dt_t = datetime.fromisoformat(t).replace(tzinfo=timezone.utc)
        cand.append((dt_t, i))
    pasados = [c for c in cand if c[0] <= ahora]
    _, idx = max(pasados) if pasados else min(cand)

    def val(k):
        arr = h.get(k) or []
        return arr[idx] if idx < len(arr) else None
    return val("wind_speed_10m"), val("wind_direction_10m"), val("wind_gusts_10m"), tiempos[idx]


def consultar_viento_ensamble(celdas):
    pesos = CFG["viento"]["modelos"]
    lote = CFG["viento"]["lote"]
    por_modelo = {}
    for modelo in pesos:
        res = {}
        for i in range(0, len(celdas), lote):
            sub = celdas[i:i + lote]
            params = {
                "latitude": ",".join(str(c[0]) for c in sub),
                "longitude": ",".join(str(c[1]) for c in sub),
                "hourly": "wind_speed_10m,wind_direction_10m,wind_gusts_10m",
                "wind_speed_unit": "kmh", "models": modelo, "timezone": "UTC",
                "past_hours": 1, "forecast_hours": 2,
            }
            j = _get_json(OPEN_METEO, params)
            if isinstance(j, dict):
                j = [j]
            for c, item in zip(sub, j):
                res[c] = _hora_actual(item) if isinstance(item, dict) and "hourly" in item else None
        por_modelo[modelo] = res
        log(f"  modelo {modelo}: {sum(1 for v in res.values() if v)}/{len(celdas)} celdas con dato")

    salida = {}
    for c in celdas:
        obs, hora = [], None
        for modelo, w in pesos.items():
            r = por_modelo[modelo].get(c)
            if r:
                obs.append((w, r[0], r[1], r[2]))
                hora = r[3]
        comb = combinar_modelos(obs)
        if comb:
            comb["hora_utc"] = hora
        salida[clave_celda(*c)] = comb
    return salida


def obtener_viento(celdas):
    """Devuelve (dict clave_celda->viento, estado_fuente, edad_min)."""
    if MOCK:
        mk = {"vel_kmh": 24.0, "rafaga_kmh": 38.0, "dir_desde": 270.0, "dir_hacia": 90.0,
              "disp_dir": 12.0, "disp_vel": 5.0, "n_modelos": 3,
              "hora_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:00")}
        return {clave_celda(*c): dict(mk) for c in celdas}, "MOCK", 0
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    archivo = CACHE_DIR / "viento_cache.json"
    ahora = datetime.now(timezone.utc)
    claves = [clave_celda(*c) for c in celdas]
    cache = None
    if archivo.exists():
        try:
            cache = json.loads(archivo.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            cache = None

    def edad(c):
        return (ahora - datetime.fromisoformat(c["generado_utc"])).total_seconds() / 60

    ttl = CFG["viento"]["ttl_minutos"]
    if cache and edad(cache) < ttl and all(k in cache["celdas"] for k in claves):
        log(f"Viento desde caché (edad {edad(cache):.0f} min, TTL {ttl} min)")
        return cache["celdas"], "CACHE", edad(cache)
    try:
        log(f"Consultando Open-Meteo para {len(celdas)} celdas...")
        nuevo = consultar_viento_ensamble(celdas)
        if sum(1 for v in nuevo.values() if v) == 0:
            raise RuntimeError("Open-Meteo no devolvió datos de viento.")
        archivo.write_text(json.dumps({"generado_utc": ahora.isoformat(), "celdas": nuevo}), encoding="utf-8")
        return nuevo, "OK", 0
    except Exception as e:  # noqa: BLE001
        advertir(f"Falló la consulta de viento: {e}")
        if cache and edad(cache) < CFG["viento"]["max_antiguedad_cache_min"]:
            advertir(f"Se usa el último viento guardado (edad {edad(cache):.0f} min).")
            return cache["celdas"], "DESACTUALIZADO", edad(cache)
        return {}, "SIN_DATO", None


# ---------------------------------------------------------------------------
# Modo prueba (solo con MOCK=1): datos sintéticos para ver la capa sin red
# ---------------------------------------------------------------------------
def mock_incendios(est):
    def desplazar(e, rumbo, km):
        dlat = km * math.cos(math.radians(rumbo)) / 111.2
        dlon = km * math.sin(math.radians(rumbo)) / (111.2 * math.cos(math.radians(e["lat"])))
        return e["lat"] + dlat, e["lon"] + dlon
    pruebas = [(est[0], 270, 7, "PRUEBA Oeste 7 km"), (est[len(est) // 2], 45, 18, "PRUEBA NE 18 km"),
               (est[-1], 180, 4, "PRUEBA Sur 4 km")]
    out = []
    for e, rumbo, km, nom in pruebas:
        la, lo = desplazar(e, rumbo, km)
        out.append({"nombre": nom, "estado": "En combate", "superficie_ha": 120, "fecha": None, "lat": la, "lon": lo})
    return out


# ---------------------------------------------------------------------------
# Cruce y salidas
# ---------------------------------------------------------------------------
def evaluar(est, fuegos, fuente_incendios_ok, viento):
    crit = CFG["criterios"]
    ahora = datetime.now(timezone.utc)
    stamp_utc = ahora.isoformat(timespec="minutes")
    stamp_cl = ahora.astimezone(TZ_CL).strftime("%Y-%m-%d %H:%M")
    puntos, lineas = [], []
    for e in est:
        w = viento.get(clave_celda(*celda_de(e["lat"], e["lon"])))
        props = {
            "id": e["id"], "nombre": e["nombre"], "tipo": e["tipo"], "comuna": e["comuna"],
            "viento_kmh": None, "rafaga_kmh": None, "dir_desde_deg": None, "dir_hacia_deg": None,
            "dir_desde_card": "", "disp_dir_deg": None, "viento_hora_utc": None,
            "estado": "SIN_INCENDIO", "tipo_viento": "", "fuego_nombre": "", "fuego_estado": "",
            "fuego_dist_km": None, "fuego_rumbo_deg": None, "delta_deg": None, "n_incendios_radio": 0,
            "fuente_incendios_ok": fuente_incendios_ok, "actualizado_utc": stamp_utc, "actualizado_cl": stamp_cl,
        }
        if w:
            props.update({
                "viento_kmh": round(w["vel_kmh"], 1),
                "rafaga_kmh": None if w.get("rafaga_kmh") is None else round(w["rafaga_kmh"], 1),
                "dir_desde_deg": None if w["dir_desde"] is None else round(w["dir_desde"]),
                "dir_hacia_deg": None if w["dir_hacia"] is None else round(w["dir_hacia"]),
                "dir_desde_card": cardinal(w["dir_desde"]),
                "disp_dir_deg": None if w.get("disp_dir") is None else round(w["disp_dir"]),
                "viento_hora_utc": w.get("hora_utc"),
            })
        mejor = None
        cercanos = 0
        for f in fuegos:
            d = haversine_km(e["lat"], e["lon"], f["lat"], f["lon"])
            if d > crit["radio_busqueda_km"]:
                continue
            cercanos += 1
            rumbo = rumbo_deg(e["lat"], e["lon"], f["lat"], f["lon"])   # del establecimiento hacia el fuego
            delta = dif_ang(w["dir_desde"], rumbo) if (w and w["dir_desde"] is not None) else None
            est_f = clasificar(d, delta, w["vel_kmh"] if w else None, crit) if w else "SIN_DATO_VIENTO"
            clave = (RANK.get(est_f, 1.5), -d)
            if mejor is None or clave > mejor[0]:
                mejor = (clave, f, d, rumbo, delta, est_f)
        props["n_incendios_radio"] = cercanos
        if mejor:
            _, f, d, rumbo, delta, est_f = mejor
            props.update({
                "estado": est_f, "tipo_viento": tipo_viento(delta, crit["tolerancia_grados"]),
                "fuego_nombre": f["nombre"], "fuego_estado": f["estado"], "fuego_dist_km": round(d, 1),
                "fuego_rumbo_deg": round(rumbo), "delta_deg": None if delta is None else round(delta),
            })
            lineas.append({"type": "Feature",
                           "geometry": {"type": "LineString", "coordinates": [[round(e["lon"], 5), round(e["lat"], 5)],
                                                                                [round(f["lon"], 5), round(f["lat"], 5)]]},
                           "properties": {"id": e["id"], "nombre": e["nombre"], "estado": est_f,
                                          "fuego_nombre": f["nombre"], "dist_km": round(d, 1),
                                          "tipo_viento": props["tipo_viento"], "actualizado_cl": stamp_cl}})
        elif not fuente_incendios_ok:
            props["estado"] = "SIN_DATO_CONAF"
        puntos.append({"type": "Feature",
                       "geometry": {"type": "Point", "coordinates": [round(e["lon"], 5), round(e["lat"], 5)]},
                       "properties": props})
    return puntos, lineas, stamp_utc, stamp_cl


def escribir_salidas(puntos, lineas, fuegos, fuente_ok, estado_viento, edad_viento, stamp_utc, stamp_cl):
    SALIDA.mkdir(parents=True, exist_ok=True)

    def gj(feats):
        return json.dumps({"type": "FeatureCollection", "features": feats}, ensure_ascii=False)
    (SALIDA / "establecimientos_viento.geojson").write_text(gj(puntos), encoding="utf-8")
    (SALIDA / "lineas_fuego_establecimiento.geojson").write_text(gj(lineas), encoding="utf-8")
    conteo = {}
    for p in puntos:
        k = p["properties"]["estado"]
        conteo[k] = conteo.get(k, 0) + 1
    estado = {
        "actualizado_utc": stamp_utc, "actualizado_cl": stamp_cl,
        "establecimientos": len(puntos), "incendios_activos_en_area": len(fuegos),
        "fuente_incendios_ok": fuente_ok, "fuente_viento": estado_viento,
        "edad_viento_min": None if edad_viento is None else round(edad_viento),
        "conteo_estados": conteo, "advertencias": ADVERTENCIAS,
    }
    (SALIDA / "estado.json").write_text(json.dumps(estado, ensure_ascii=False, indent=2), encoding="utf-8")
    filas = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in sorted(conteo.items(), key=lambda kv: -RANK.get(kv[0], 0)))
    adv = "".join(f"<li>{a}</li>" for a in ADVERTENCIAS) or "<li>Sin advertencias</li>"
    html = f"""<!doctype html><html lang="es"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Monitor viento–incendios</title>
<style>body{{font-family:system-ui,sans-serif;max-width:640px;margin:24px auto;padding:0 16px;color:#1B2140}}
td,th{{padding:4px 12px;border-bottom:1px solid #ddd;text-align:left}}</style></head><body>
<h2>Monitor viento–incendios · Servicio de Salud Maule</h2>
<p>Última actualización: <b>{stamp_cl}</b> (hora Chile). Fuente CONAF: <b>{'OK' if fuente_ok else 'SIN DATO'}</b>.
Viento: <b>{estado_viento}</b>. Incendios activos en el área: <b>{len(fuegos)}</b>.</p>
<table><tr><th>Estado</th><th>Establecimientos</th></tr>{filas}</table>
<h3>Capas</h3><ul>
<li><a href="establecimientos_viento.geojson">establecimientos_viento.geojson</a></li>
<li><a href="lineas_fuego_establecimiento.geojson">lineas_fuego_establecimiento.geojson</a></li>
<li><a href="estado.json">estado.json</a></li></ul>
<h3>Advertencias</h3><ul>{adv}</ul>
<p style="color:#5B6280;font-size:13px">Herramienta de tamizaje basada en modelos de 9–25 km. No reemplaza el análisis de CONAF/SENAPRED.</p>
</body></html>"""
    (SALIDA / "index.html").write_text(html, encoding="utf-8")
    log(f"Salidas escritas en {SALIDA}: {conteo}")


def main():
    est = cargar_establecimientos()
    resultado = cargar_incendios(est)
    fuegos, fuente_ok = resultado
    celdas = sorted({celda_de(e["lat"], e["lon"]) for e in est})
    log(f"Celdas únicas de viento: {len(celdas)} (de {len(est)} establecimientos)")
    viento, estado_viento, edad = obtener_viento(celdas)
    puntos, lineas, s_utc, s_cl = evaluar(est, fuegos, fuente_ok, viento)
    escribir_salidas(puntos, lineas, fuegos, fuente_ok, estado_viento, edad, s_utc, s_cl)


if __name__ == "__main__":
    main()
