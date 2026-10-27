#!/usr/bin/env python3
"""
Build the venue layout library (schema vsp-layout/1) for every stadium in venues.csv.

For each venue it pulls OpenStreetMap data from the Overpass API and writes
venues/<slug>.geojson with:
  - mapped features (stadium outline, playing field, gates, elevators, nearby
    enclosed buildings as shelter candidates) -- all marked unverified
  - modeled seating zones and concourse segments generated from the gap between
    the field and the stadium's outer outline -- marked modeled=yes

Hand-traced features in venues/edits/<slug>.geojson are merged in on every run and
replace the modeled zones/concourses of the same kind, so tracing work is never lost.

Also writes venues/index.json (for the dashboard's stadium picker) and
venues/coverage.csv (what was found for each venue).

Standard library only.  Usage:
  python build_venues.py                       # all venues
  python build_venues.py --only lane-stadium   # one or more slugs (comma list)
  python build_venues.py --limit 5             # first N (testing)
  python build_venues.py --cache cache_dir     # read/write raw Overpass responses
"""
import argparse, csv, json, math, os, re, sys, time, unicodedata, urllib.parse, urllib.request
from datetime import datetime, timezone

SCHEMA = 'vsp-layout/1'
MODEL_ID = 'vsp-model/1'
UA = 'stadium-venue-library builder (github.com/zeringueil-source/stadium-venue-library)'
OVERPASS = ['https://overpass-api.de/api/interpreter', 'https://overpass.kumi.systems/api/interpreter']
R_EARTH = 6371008.8
M2_PER_FT2 = 0.09290304

# Model parameters (documented in the README; edit here)
APRON_M = 8.0              # gap between field edge and first seats
CONCOURSE_MAX_M = 12.0     # depth of the modeled concourse band at the outer edge
CONCOURSE_FRAC = 0.25      # ...but never more than this share of the field-to-wall depth
MIN_BOWL_DEPTH_M = 12.0    # rays with less field-to-wall depth than this get no seats
EGRESS_TARGET_MIN = 8.0    # modeled exit width assumes each zone clears in this time
STEPPED_RATE = 20.0        # persons / ft / min on stepped surfaces (Green Guide 66 / m / min)
SECTORS = 8
RAY_STEP_DEG = 2
BUILDING_RADIUS_M = 350
BUILDING_MIN_M2 = 400
MAX_BUILDINGS = 30
SKIP_BUILDINGS = {'roof', 'carport', 'shed', 'hut', 'tent', 'ruins', 'construction', 'kiosk', 'grandstand',
                  'stadium', 'stands', 'container', 'shelter', 'toilets', 'service', 'no',
                  # private homes are not public shelters
                  'house', 'apartments', 'residential', 'detached', 'semidetached_house', 'terrace', 'bungalow',
                  'dormitory', 'garage', 'garages', 'farm', 'cabin', 'static_caravan', 'houseboat'}
COMPASS = ['North', 'Northeast', 'East', 'Southeast', 'South', 'Southwest', 'West', 'Northwest']


# ----------------------------------------------------------------------------- utils
def slugify(s):
    s = ''.join(ch for ch in unicodedata.normalize('NFKD', str(s)) if not unicodedata.combining(ch)).lower()
    return re.sub(r'[^a-z0-9]+', '-', s).strip('-') or 'venue'


def http(url, data=None, timeout=180):
    req = urllib.request.Request(url, data=data.encode() if data else None, headers={'User-Agent': UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8')


class Proj:
    """Local equirectangular projection in metres around (lat0, lon0)."""
    def __init__(self, lat0, lon0):
        self.lat0, self.lon0 = lat0, lon0
        self.kx = math.cos(math.radians(lat0)) * math.pi / 180 * R_EARTH
        self.ky = math.pi / 180 * R_EARTH

    def fwd(self, lon, lat):
        return ((lon - self.lon0) * self.kx, (lat - self.lat0) * self.ky)

    def inv(self, x, y):
        return (round(self.lon0 + x / self.kx, 7), round(self.lat0 + y / self.ky, 7))


def ring_area(pts):  # signed, metres^2 (projected points)
    return sum(pts[i][0] * pts[(i + 1) % len(pts)][1] - pts[(i + 1) % len(pts)][0] * pts[i][1] for i in range(len(pts))) / 2


def ring_centroid(pts):
    a = ring_area(pts)
    if abs(a) < 1e-9:
        return (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))
    cx = cy = 0.0
    for i in range(len(pts)):
        x0, y0 = pts[i]; x1, y1 = pts[(i + 1) % len(pts)]
        f = x0 * y1 - x1 * y0
        cx += (x0 + x1) * f; cy += (y0 + y1) * f
    return (cx / (6 * a), cy / (6 * a))


def point_in(pt, pts):
    x, y = pt; inside = False; j = len(pts) - 1
    for i in range(len(pts)):
        xi, yi = pts[i]; xj, yj = pts[j]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / ((yj - yi) or 1e-12) + xi:
            inside = not inside
        j = i
    return inside


def dist_to_ring(pt, pts):
    if point_in(pt, pts):
        return 0.0
    best = float('inf')
    for i in range(len(pts)):
        ax, ay = pts[i]; bx, by = pts[(i + 1) % len(pts)]
        dx, dy = bx - ax, by - ay; L2 = dx * dx + dy * dy
        t = 0 if L2 == 0 else max(0, min(1, ((pt[0] - ax) * dx + (pt[1] - ay) * dy) / L2))
        best = min(best, math.hypot(ax + t * dx - pt[0], ay + t * dy - pt[1]))
    return best


def ray_max_hit(c, ang, pts):
    """Farthest intersection distance of a ray from c at math-angle ang with polygon pts (None if none)."""
    ux, uy = math.cos(ang), math.sin(ang); best = None
    for i in range(len(pts)):
        ax, ay = pts[i]; bx, by = pts[(i + 1) % len(pts)]
        ex, ey = bx - ax, by - ay
        den = ux * ey - uy * ex
        if abs(den) < 1e-12:
            continue
        wx, wy = ax - c[0], ay - c[1]
        t = (wx * ey - wy * ex) / den     # along ray
        s = (wx * uy - wy * ux) / den     # along edge
        if t > 0 and 0 <= s <= 1 and (best is None or t > best):
            best = t
    return best


def principal_axis(pts):
    """Angle (radians, math convention) of the long axis of a point set."""
    n = len(pts); mx = sum(p[0] for p in pts) / n; my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts); syy = sum((p[1] - my) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    return 0.5 * math.atan2(2 * sxy, sxx - syy)


# ----------------------------------------------------------------------------- OSM
def overpass_query(lat, lon):
    r, rb = 700, BUILDING_RADIUS_M + 300
    return f"""[out:json][timeout:90];
(
  way(around:{r},{lat},{lon})["leisure"="stadium"];
  relation(around:{r},{lat},{lon})["leisure"="stadium"];
  way(around:{r},{lat},{lon})["leisure"="pitch"];
  relation(around:{r},{lat},{lon})["leisure"="pitch"];
  way(around:{rb},{lat},{lon})["building"];
  relation(around:{rb},{lat},{lon})["building"];
  node(around:{r},{lat},{lon})["entrance"];
  node(around:{r},{lat},{lon})["highway"="elevator"];
);
out geom;"""


def fetch_osm(lat, lon, cache_path=None):
    if cache_path and os.path.exists(cache_path):
        with open(cache_path) as f:
            return json.load(f), True
    q = overpass_query(lat, lon); last = None
    for attempt in range(4):
        ep = OVERPASS[attempt % len(OVERPASS)]
        try:
            js = json.loads(http(ep, 'data=' + urllib.parse.quote(q)))
            if cache_path:
                with open(cache_path, 'w') as f:
                    json.dump(js, f)
            return js, False
        except Exception as e:  # 429 / 504 / timeouts: back off and retry
            last = e
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f'Overpass failed: {last}')


def geocode(name, location):
    q = urllib.parse.urlencode({'q': f'{name}, {location}', 'format': 'jsonv2', 'limit': 1, 'countrycodes': 'us'})
    js = json.loads(http('https://nominatim.openstreetmap.org/search?' + q))
    time.sleep(1.1)  # Nominatim policy: max 1 request per second
    return (float(js[0]['lat']), float(js[0]['lon'])) if js else None


def assemble_rings(members):
    """Join relation member ways (role outer) into closed rings of (lon, lat)."""
    segs = [[(g['lon'], g['lat']) for g in m['geometry']] for m in members
            if m.get('type') == 'way' and m.get('role', 'outer') in ('outer', '') and m.get('geometry')]
    rings = []
    while segs:
        ring = segs.pop(0)
        changed = True
        while ring[0] != ring[-1] and changed:
            changed = False
            for i, s in enumerate(segs):
                if s[0] == ring[-1]: ring += s[1:]
                elif s[-1] == ring[-1]: ring += s[::-1][1:]
                elif s[-1] == ring[0]: ring = s[:-1] + ring
                elif s[0] == ring[0]: ring = s[::-1][:-1] + ring
                else: continue
                segs.pop(i); changed = True; break
        if ring[0] == ring[-1] and len(ring) >= 4:
            rings.append(ring)
    return rings


def osm_polygons(js):
    """Yield (element, outer_ring_lonlat) for closed ways and multipolygon relations."""
    for el in js.get('elements', []):
        if el['type'] == 'way' and el.get('geometry'):
            c = [(g['lon'], g['lat']) for g in el['geometry']]
            if len(c) >= 4 and c[0] == c[-1]:
                yield el, c
        elif el['type'] == 'relation' and el.get('members'):
            rings = assemble_rings(el['members'])
            if rings:
                yield el, max(rings, key=lambda r: abs(ring_area(r)))


# ----------------------------------------------------------------------------- build
def feature(ftype, geom, **props):
    p = {'ftype': ftype}
    p.update({k: v for k, v in props.items() if v not in (None, '')})
    return {'type': 'Feature', 'geometry': geom, 'properties': p}


def poly_geom(ring_lonlat):
    r = [list(c) for c in ring_lonlat]
    if r[0] != r[-1]:
        r.append(r[0])
    return {'type': 'Polygon', 'coordinates': [r]}


def build_venue(v, js, run_date):
    lat, lon = v['lat'], v['lon']
    P = Proj(lat, lon)
    venue_pt = (0.0, 0.0)
    feats, notes = [], []
    osm_src = lambda el: f"OpenStreetMap {el['type']}/{el['id']} (ODbL)"

    polys = []
    for el, ring in osm_polygons(js):
        pts = [P.fwd(*c) for c in ring[:-1]]
        if len(pts) < 3:
            continue
        polys.append({'el': el, 'tags': el.get('tags', {}), 'ring': ring, 'pts': pts, 'area': abs(ring_area(pts)),
                      'dist': dist_to_ring(venue_pt, pts), 'cen': ring_centroid(pts)})

    # --- playing field: the main pitch at the venue point (game sports preferred over courts)
    GAME = ('american_football', 'soccer', 'baseball', 'multi', 'rugby', 'football', 'lacrosse')
    pitches = [c for c in polys if c['tags'].get('leisure') == 'pitch' and c['area'] > 1500 and c['dist'] <= 80]
    game = [c for c in pitches if any(s in (c['tags'].get('sport') or '') for s in GAME)] or pitches
    containing = [c for c in game if c['dist'] == 0]
    field = max(containing, key=lambda c: c['area']) if containing else (min(game, key=lambda c: c['dist']) if game else None)

    # --- stadium outline: the smallest stadium polygon (building=stadium or leisure=stadium) that contains
    #     the field; without a field, the one containing / nearest the venue point
    stad = [c for c in polys if (c['tags'].get('building') == 'stadium' or c['tags'].get('leisure') == 'stadium')
            and c['area'] > 4000]
    outline = None
    if field:
        hold = [c for c in stad if c is not field and point_in(field['cen'], c['pts']) and c['area'] > 1.3 * field['area']]
        outline = min(hold, key=lambda c: c['area']) if hold else None
    if not outline:
        near = [c for c in stad if c['dist'] <= 150]
        inside = [c for c in near if c['dist'] == 0]
        outline = min(inside, key=lambda c: c['area']) if inside else (min(near, key=lambda c: c['dist']) if near else None)
    shown = [c for c in stad if c['dist'] <= 150 and (c is outline or not field or point_in(field['cen'], c['pts'])
                                                       or c['dist'] <= 20)]
    for c in shown:
        feats.append(feature('structure', poly_geom(c['ring']), name=c['tags'].get('name') or 'Stadium outline',
                             verified='no', modeled='no', source=osm_src(c['el']), src_date=run_date))
    if field:
        feats.append(feature('structure', poly_geom(field['ring']), name='Playing field', verified='no', modeled='no',
                             source=osm_src(field['el']), src_date=run_date))

    # --- grandstand / stands structures
    for c in polys:
        if c['tags'].get('building') in ('grandstand', 'stands') and c['dist'] <= 250:
            feats.append(feature('structure', poly_geom(c['ring']), name=c['tags'].get('name') or 'Grandstand',
                                 verified='no', modeled='no', source=osm_src(c['el']), src_date=run_date))

    # --- modeled seating bowl + concourse
    modeled = {'zones': 0, 'concourses': 0}
    if outline:
        if field:
            fpts = field['pts']
        else:  # no mapped pitch: assume an inner field at 45% scale of the outline
            oc = outline['cen']
            fpts = [(oc[0] + (x - oc[0]) * 0.45, oc[1] + (y - oc[1]) * 0.45) for x, y in outline['pts']]
            notes.append('no playing field mapped; field assumed at 45% of outline')
        c0 = ring_centroid(fpts)
        if not point_in(c0, outline['pts']):
            notes.append('field centre outside stadium outline; seating not modeled')
        else:
            axis = principal_axis(fpts)
            half = math.pi / SECTORS  # sector boundaries offset so each sideline/end is one whole sector
            sectors = []
            for s in range(SECTORS):
                a0 = axis - half + s * 2 * half
                seat_in, seat_out, conc_in, conc_out = [], [], [], []
                steps = max(2, int(round(math.degrees(2 * half) / RAY_STEP_DEG)))
                ok = 0
                for k in range(steps + 1):
                    a = a0 + 2 * half * k / steps
                    rf = ray_max_hit(c0, a, fpts) or 0.0
                    ro = ray_max_hit(c0, a, outline['pts'])
                    if ro is None:
                        ro = rf
                    depth = ro - rf
                    if depth < MIN_BOWL_DEPTH_M:
                        r_in = r_so = r_ci = r_co = rf if rf < ro else ro
                    else:
                        cd = min(CONCOURSE_MAX_M, CONCOURSE_FRAC * depth)
                        r_in, r_so = rf + min(APRON_M, depth * 0.2), ro - cd
                        r_ci, r_co = ro - cd, ro
                        ok += 1
                    pt = lambda r: (c0[0] + r * math.cos(a), c0[1] + r * math.sin(a))
                    seat_in.append(pt(r_in)); seat_out.append(pt(r_so))
                    conc_in.append(pt(r_ci)); conc_out.append(pt(r_co))
                seat = seat_in + seat_out[::-1]
                conc = conc_in + conc_out[::-1]
                bmid = (90 - math.degrees(a0 + half)) % 360  # compass bearing of sector centre
                sectors.append({'seat': seat, 'conc': conc, 'seat_area': abs(ring_area(seat)) if ok else 0,
                                'conc_area': abs(ring_area(conc)) if ok else 0, 'name': COMPASS[int(((bmid + 22.5) % 360) // 45)]})
            total = sum(s['seat_area'] for s in sectors)
            cap = v['capacity'] or 0
            roof = (v['roof'] or '').lower()
            enclosed = roof in ('dome', 'fixed')
            used = {}
            for s in sectors:
                if s['seat_area'] < 50 or total <= 0:
                    continue
                nm = s['name']; used[nm] = used.get(nm, 0) + 1
                nm = nm if used[nm] == 1 else f'{nm} {used[nm]}'
                persons = round(cap * s['seat_area'] / total) if cap else None
                exit_ft = round(persons / (STEPPED_RATE * EGRESS_TARGET_MIN), 1) if persons else None
                ring = [P.inv(x, y) for x, y in s['seat']]
                feats.append(feature('seating_zone', poly_geom(ring), name=f'{nm} stands (modeled)', persons=persons,
                                     exit_ft=exit_ft, exit_type='stepped', verified='no', modeled='yes',
                                     source=f'{MODEL_ID}: modeled from OSM outline; capacity spread by area',
                                     src_date=run_date,
                                     notes=f'Modeled, not mapped. Exit width assumes an {int(EGRESS_TARGET_MIN)}-minute egress design. Replace with traced sections.'))
                modeled['zones'] += 1
                if s['conc_area'] > 50:
                    ring = [P.inv(x, y) for x, y in s['conc']]
                    feats.append(feature('shelter', poly_geom(ring), name=f'{nm} concourse (modeled)',
                                         shel_class='enclosed_concourse' if enclosed else 'covered_concourse',
                                         verified='no', modeled='yes', source=f'{MODEL_ID}: outer band of stadium outline',
                                         src_date=run_date,
                                         notes='Modeled, not mapped. Location and size of the real concourse are unknown; confirm with the venue.'))
                    modeled['concourses'] += 1
            if enclosed:
                notes.append('domed venue: the bowl itself is enclosed; lightning shelter is mainly for plazas, lots and tailgates')
            elif roof in ('retractable', 'skylight'):
                notes.append(f'{roof} roof: modeled plan applies when the roof is open')
    else:
        notes.append('no stadium outline in OSM; seating not modeled')

    # --- nearby enclosed buildings as shelter candidates
    ref_pts = outline['pts'] if outline else None
    cands = []
    for c in polys:
        b = c['tags'].get('building')
        if not b or b in SKIP_BUILDINGS or c['area'] < BUILDING_MIN_M2:
            continue
        if outline and (c is outline or point_in(c['cen'], outline['pts'])):
            continue
        d = min(dist_to_ring(p, ref_pts) for p in c['pts']) if ref_pts else c['dist']
        if d <= BUILDING_RADIUS_M:
            cands.append((d, c))
    cands.sort(key=lambda t: t[0])
    for d, c in cands[:MAX_BUILDINGS]:
        t = c['tags']; b = t.get('building')
        note = 'Auto-drafted from OSM. Confirm it is open, unlocked and usable during events.'
        if b in ('parking',) or t.get('amenity') == 'parking':
            note = 'Parking structure. Confirm with the venue whether it is an approved lightning shelter.'
        feats.append(feature('shelter', poly_geom(c['ring']), name=t.get('name') or (b.replace('_', ' ').title() if b != 'yes' else 'Building'),
                             shel_class='enclosed_building', verified='no', modeled='no', source=osm_src(c['el']),
                             src_date=run_date, notes=note, level=t.get('building:levels')))

    # --- gates and elevators near the stadium
    n_gates = n_elev = 0
    for el in js.get('elements', []):
        if el['type'] != 'node':
            continue
        t = el.get('tags', {}); pt = P.fwd(el['lon'], el['lat'])
        near = (dist_to_ring(pt, outline['pts']) <= 120) if outline else (math.hypot(*pt) <= 400)
        if not near:
            continue
        g = {'type': 'Point', 'coordinates': [el['lon'], el['lat']]}
        if t.get('highway') == 'elevator':
            feats.append(feature('vertical', g, name=t.get('name') or 'Elevator', vkind='elevator', ada='yes',
                                 verified='no', modeled='no', source=osm_src(el), src_date=run_date)); n_elev += 1
        elif t.get('entrance'):
            feats.append(feature('gate', g, name=t.get('name') or t.get('ref') or ('Gate' if t['entrance'] in ('main', 'yes') else t['entrance'].title()),
                                 ada='yes' if t.get('wheelchair') == 'yes' else None, verified='no', modeled='no',
                                 source=osm_src(el), src_date=run_date)); n_gates += 1

    cov = {'outline': bool(outline), 'outline_kind': ('building' if outline['tags'].get('building') == 'stadium' else 'leisure') if outline else '',
           'field': bool(field), 'modeled_zones': modeled['zones'], 'modeled_concourses': modeled['concourses'],
           'buildings': min(len(cands), MAX_BUILDINGS), 'gates': n_gates, 'elevators': n_elev, 'notes': '; '.join(notes)}
    return feats, cov


def assign_ids(feats):
    """Stable ids so hand edits can be matched across runs: osm-<type>-<id> or model-<name>."""
    seen = {}
    for f in feats:
        p = f['properties']
        m = re.match(r'OpenStreetMap (\w+)/(\d+)', p.get('source', ''))
        base = f'osm-{m.group(1)}-{m.group(2)}' if m else ('model-' + slugify(p.get('name', 'feature')) if p.get('modeled') == 'yes' else 'gen-' + slugify(p.get('name', 'feature')))
        if p['ftype'] == 'structure' and p.get('name') == 'Playing field':
            base += '-field'
        seen[base] = seen.get(base, 0) + 1
        f['id'] = base if seen[base] == 1 else f'{base}-{seen[base]}'
    return feats


def merge_edits(feats, edits_path):
    """Merge venues/edits/<slug>.geojson (exported from the dashboard or traced in GIS).

    - an edited feature replaces the generated feature with the same id
    - generated ids listed in the edits file's "removed_ids" are dropped
    - if the edits contain traced (non-modeled, non-OSM) seating zones, the modeled seating zones are dropped;
      same for traced concourse shelters vs modeled concourses
    """
    if not os.path.exists(edits_path):
        return feats, 0
    with open(edits_path, encoding='utf-8') as f:
        ed = json.load(f)
    efeats = [x for x in ed.get('features', []) if (x.get('properties') or {}).get('source') != 'example' and x.get('geometry')]
    removed = set(ed.get('removed_ids') or [])
    eids = {str(x.get('id')) for x in efeats if x.get('id') is not None}
    yes = lambda v: str(v).lower() in ('yes', 'true', '1')
    traced = [x for x in efeats if not yes((x.get('properties') or {}).get('modeled'))
              and not str((x.get('properties') or {}).get('source', '')).startswith('OpenStreetMap')]
    drop_zones = any((x['properties'] or {}).get('ftype') == 'seating_zone' for x in traced)
    drop_conc = any((x['properties'] or {}).get('ftype') == 'shelter' and 'concourse' in str((x['properties'] or {}).get('shel_class', '')) for x in traced)
    keep = []
    for x in feats:
        p = x['properties']
        if x['id'] in eids or x['id'] in removed:
            continue
        if p.get('modeled') == 'yes' and ((p['ftype'] == 'seating_zone' and drop_zones) or (p['ftype'] == 'shelter' and drop_conc)):
            continue
        keep.append(x)
    return keep + efeats, len(efeats)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default='venues.csv')
    ap.add_argument('--out', default='venues')
    ap.add_argument('--only', default='')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--cache', default='')
    ap.add_argument('--sleep', type=float, default=4.0)
    a = ap.parse_args()

    with open(a.csv, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    names = [r['name'] for r in rows]
    for r in rows:
        base = slugify(r['name'])
        r['slug'] = base if names.count(r['name']) == 1 else slugify(r['name'] + ' ' + r['location'].split(',')[0])
        r['capacity'] = int(r['capacity']) if r.get('capacity') else None
        r['lat'] = float(r['lat']) if r.get('lat') else None
        r['lon'] = float(r['lon']) if r.get('lon') else None
    only = {s.strip() for s in a.only.split(',') if s.strip()}
    todo = [r for r in rows if not only or r['slug'] in only]
    if a.limit:
        todo = todo[:a.limit]

    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.join(a.out, 'edits'), exist_ok=True)
    if a.cache:
        os.makedirs(a.cache, exist_ok=True)
    run_date = datetime.now(timezone.utc).strftime('%Y-%m-%d')

    index_path = os.path.join(a.out, 'index.json')
    index = {}
    if os.path.exists(index_path):
        with open(index_path) as f:
            index = {v['slug']: v for v in json.load(f).get('venues', [])}

    failures = streak = ok = 0
    for i, v in enumerate(todo, 1):
        cached = True
        try:
            if v['lat'] is None:
                g = geocode(v['name'], v['location'])
                if not g:
                    raise RuntimeError('could not geocode')
                v['lat'], v['lon'] = g
            js, cached = fetch_osm(v['lat'], v['lon'], os.path.join(a.cache, v['slug'] + '.json') if a.cache else None)
            feats, cov = build_venue(v, js, run_date)
            feats = assign_ids(feats)
            feats, n_ed = merge_edits(feats, os.path.join(a.out, 'edits', v['slug'] + '.geojson'))
            cov['traced_features'] = n_ed
            fc = {'type': 'FeatureCollection', 'schema': SCHEMA,
                  'venue': {'name': v['name'], 'slug': v['slug'], 'lat': v['lat'], 'lon': v['lon'], 'capacity': v['capacity'],
                            'location': v['location'], 'leagues': v['leagues'], 'teams': v['teams'], 'roof': v['roof'] or None},
                  'generated': datetime.now(timezone.utc).isoformat(timespec='seconds'),
                  'attribution': 'Map data © OpenStreetMap contributors (ODbL). Modeled features: ' + MODEL_ID,
                  'coverage': cov, 'features': feats}
            for k, ft in enumerate(fc['features']):
                if ft.get('id') is None:
                    ft['id'] = f"edit-{k}"
            with open(os.path.join(a.out, v['slug'] + '.geojson'), 'w', encoding='utf-8') as f:
                json.dump(fc, f, ensure_ascii=False, separators=(',', ':'))
            index[v['slug']] = {'slug': v['slug'], 'name': v['name'], 'location': v['location'], 'leagues': v['leagues'],
                                'teams': v['teams'], 'capacity': v['capacity'], 'roof': v['roof'] or None,
                                'lat': round(v['lat'], 6), 'lon': round(v['lon'], 6), 'coverage': cov, 'updated': run_date}
            print(f"[{i}/{len(todo)}] {v['slug']}: outline={cov['outline']} field={cov['field']} zones={cov['modeled_zones']} "
                  f"bldgs={cov['buildings']} gates={cov['gates']} {cov['notes']}", flush=True)
            streak = 0; ok += 1
        except Exception as e:
            failures += 1; streak += 1
            print(f"[{i}/{len(todo)}] {v['slug']}: FAILED {e}", flush=True)
            if streak >= 8:
                print('8 failures in a row: Overpass looks unavailable, stopping this run early (existing files kept).')
                break
        if not cached:
            time.sleep(a.sleep)

    with open(index_path, 'w', encoding='utf-8') as f:
        json.dump({'schema': SCHEMA, 'generated': datetime.now(timezone.utc).isoformat(timespec='seconds'),
                   'count': len(index), 'venues': sorted(index.values(), key=lambda x: x['name'])}, f, ensure_ascii=False, indent=0)
    with open(os.path.join(a.out, 'coverage.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['slug', 'name', 'leagues', 'outline', 'field', 'modeled_zones', 'modeled_concourses', 'buildings', 'gates', 'elevators', 'traced_features', 'notes'])
        for x in sorted(index.values(), key=lambda x: x['name']):
            c = x['coverage']
            w.writerow([x['slug'], x['name'], x['leagues'], c['outline'], c['field'], c['modeled_zones'], c['modeled_concourses'],
                        c['buildings'], c['gates'], c['elevators'], c.get('traced_features', 0), c['notes']])
    print(f'done: {ok} ok, {failures} failed, {len(todo) - ok - failures} not attempted, index has {len(index)} venues')
    # Fail loudly only if most of the run failed (e.g. Overpass down), so partial outages don't wipe nothing but still show red
    if todo and (failures > len(todo) / 2 or streak >= 8):
        sys.exit(1)


if __name__ == '__main__':
    main()
