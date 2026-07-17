#!/usr/bin/env python3
"""
Convert ABAQUS .inp (non-assembly / flat format) to MyFEM JSON.

Element mapping:
    T3D2 -> Bar3D  (material: bar)
    CPS4 -> Q4_PS  (material: plane_stress)
    CPE4 -> Q4_PE  (material: plane_strain)
    C3D8 -> H8     (material: solid3d)

BC / loads from ABAQUS keywords:
    *Boundary   ENCASTRE / PINNED / XSYMM / YSYMM / ZSYMM / numeric dof ranges
    *Cload      per-node concentrated force
    *Surface + *Dsload (TRVEC)   surface traction  -> surface load

Original element numbering is preserved (elements are NOT renumbered),
so surface loads can reference original ABAQUS element ids.

Body force via top-of-file comment:
    **BODYFORCE gx gy gz

Usage:
    python inp2json.py <in.inp> [out.json]
"""

import re
import sys
import json
import math


ELEMENT_MAP = {
    'T3D2': ('Bar3D', 'bar',          3, 2),
    'CPS4': ('Q4_PS', 'plane_stress', 2, 4),
    'CPE4': ('Q4_PE', 'plane_strain', 2, 4),
    'C3D8': ('H8',    'solid3d',      3, 8),
}

BC_TYPES = {
    'ENCASTRE': [1, 2, 3],
    'PINNED':   [1, 2, 3],
    'XSYMM':    [1],
    'YSYMM':    [2],
    'ZSYMM':    [3],
}

DOF_NAME = {1: 'x', 2: 'y', 3: 'z'}


# --------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------
def parse_inp(path):
    data = dict(
        nodes={}, elements={}, elem_type=None,
        nsets={}, elsets={},
        material=dict(E=None, nu=None, rho=0.0),
        thk=None,
        boundary=[],
        cload=[],
        surfaces={},     # name -> list of (element_id, S_id)
        dsload=[],       # list of (surface_name, type, magnitude, [dir])
        body_force=None,
    )

    with open(path) as f:
        lines = f.readlines()

    for ln in lines:
        s = ln.strip()
        if s.lower().startswith('**bodyforce'):
            p = s.split()
            if len(p) >= 4:
                data['body_force'] = [float(p[1]), float(p[2]), float(p[3])]
            break

    def read_block(start):
        rows, j = [], start
        while j < len(lines):
            s = lines[j].strip()
            if s.startswith('*') and not s.startswith('**'):
                break
            if s and not s.startswith('**'):
                rows.append([t.strip() for t in s.split(',') if t.strip()])
            j += 1
        return rows, j

    def expand_set(rows, generate):
        ids = []
        for r in rows:
            if generate and len(r) >= 3:
                s0, e0, st = int(r[0]), int(r[1]), int(r[2])
                ids += list(range(s0, e0 + 1, st))
            else:
                ids += [int(x) for x in r]
        return ids

    i = 0
    while i < len(lines):
        raw = lines[i].strip()
        if not raw or raw.startswith('**'):
            i += 1
            continue
        low = raw.lower()

        if low.startswith('*node'):
            rows, i = read_block(i + 1)
            for r in rows:
                data['nodes'][int(r[0])] = [float(x) for x in r[1:]]
            continue

        if low.startswith('*element'):
            m = re.search(r'type\s*=\s*(\w+)', raw, re.IGNORECASE)
            if m:
                data['elem_type'] = m.group(1).upper()
            rows, i = read_block(i + 1)
            for r in rows:
                data['elements'][int(r[0])] = [int(x) for x in r[1:]]
            continue

        if low.startswith('*nset'):
            m = re.search(r'nset\s*=\s*([^,\s]+)', raw, re.IGNORECASE)
            if m:
                name = m.group(1).lower()
                gen = 'generate' in low
                rows, i = read_block(i + 1)
                data['nsets'].setdefault(name, [])
                data['nsets'][name] += expand_set(rows, gen)
            else:
                i += 1
            continue

        if low.startswith('*elset'):
            m = re.search(r'elset\s*=\s*([^,\s]+)', raw, re.IGNORECASE)
            if m:
                name = m.group(1).lower()
                gen = 'generate' in low
                rows, i = read_block(i + 1)
                data['elsets'].setdefault(name, [])
                data['elsets'][name] += expand_set(rows, gen)
            else:
                i += 1
            continue

        if low.startswith('*surface'):
            m = re.search(r'name\s*=\s*([^,\s]+)', raw, re.IGNORECASE)
            name = m.group(1).lower() if m else None
            rows, i = read_block(i + 1)
            if name is not None:
                data['surfaces'].setdefault(name, [])
                # each row: <elset_or_elem>, S<k>
                data['surfaces'][name].append(('RAW', rows))
            continue

        if low.startswith('*material'):
            i += 1
            continue

        if low.startswith('*elastic'):
            rows, i = read_block(i + 1)
            if rows:
                data['material']['E'] = float(rows[0][0])
                if len(rows[0]) > 1:
                    data['material']['nu'] = float(rows[0][1])
            continue

        if low.startswith('*density'):
            rows, i = read_block(i + 1)
            if rows:
                data['material']['rho'] = float(rows[0][0])
            continue

        if low.startswith('*solid section'):
            m = re.search(r'thk\s*=\s*([\d.eE+-]+)', raw, re.IGNORECASE)
            if m:
                data['thk'] = float(m.group(1))
            rows, i = read_block(i + 1)
            if data['thk'] is None and rows:
                try:
                    data['thk'] = float(rows[0][0])
                except (ValueError, IndexError):
                    pass
            continue

        if low.startswith('*boundary'):
            rows, i = read_block(i + 1)
            data['boundary'] += rows
            continue

        if low.startswith('*cload'):
            rows, i = read_block(i + 1)
            for r in rows:
                if len(r) >= 3:
                    data['cload'].append((r[0].lower(), int(r[1]), float(r[2])))
            continue

        if low.startswith('*dsload') or low.startswith('*dload'):
            rows, i = read_block(i + 1)
            for r in rows:
                # TRVEC:  surf, TRVEC, mag, n1, n2, n3
                # P:      surf, P, mag
                if len(r) < 2:
                    continue
                surf = r[0].lower()
                ltype = r[1].upper()
                if ltype == 'TRVEC' and len(r) >= 6:
                    mag = float(r[2])
                    direction = [float(r[3]), float(r[4]), float(r[5])]
                    data['dsload'].append((surf, 'TRVEC', mag, direction))
                elif ltype.startswith('P'):
                    mag = float(r[2]) if len(r) > 2 else 0.0
                    data['dsload'].append((surf, 'P', mag, None))
                else:
                    print("  [warn] unsupported *Dsload type: %s" % ltype)
            continue

        i += 1

    return data


# --------------------------------------------------------------------
# Resolve surfaces -> list of (element_id, face_number_1based)
# --------------------------------------------------------------------
def resolve_surface_faces(data, surf_name):
    """Return list of (element_id, face_1based) for a named surface."""
    faces = []
    entries = data['surfaces'].get(surf_name, [])
    for kind, rows in entries:
        for r in rows:
            if len(r) < 2:
                continue
            target = r[0].lower()
            s_id = r[1].upper()          # like 'S3'
            m = re.match(r'S(\d+)', s_id)
            if not m:
                print("  [warn] bad face id '%s' in surface %s" % (s_id, surf_name))
                continue
            face_1based = int(m.group(1))  # ABAQUS Sn -> our 1-based face n

            # target may be an elset name or a single element id
            if target in data['elsets']:
                elem_ids = data['elsets'][target]
            else:
                try:
                    elem_ids = [int(r[0])]
                except ValueError:
                    print("  [warn] surface %s target '%s' not found"
                          % (surf_name, r[0]))
                    continue
            for eid in elem_ids:
                faces.append((eid, face_1based))
    return faces


# --------------------------------------------------------------------
# Build JSON model
# --------------------------------------------------------------------
def build_json(data):
    aba = data['elem_type']
    if aba not in ELEMENT_MAP:
        raise ValueError("Unsupported ABAQUS element type: %s" % aba)
    elem_type, mat_type, dim, nnode = ELEMENT_MAP[aba]

    # --- nodes (require contiguous 1..N, matching the C++ reader) ---
    node_ids = sorted(data['nodes'].keys())
    if node_ids != list(range(1, len(node_ids) + 1)):
        raise ValueError(
            "Node ids must be contiguous 1..N (got %d nodes, min=%d max=%d)."
            % (len(node_ids), node_ids[0], node_ids[-1]))

    nodes_json = []
    for nid in node_ids:
        c = data['nodes'][nid]
        entry = {"id": nid, "x": c[0], "y": c[1] if len(c) > 1 else 0.0}
        if dim == 3:
            entry["z"] = c[2] if len(c) > 2 else 0.0
        nodes_json.append(entry)

    # --- material ---
    mat = data['material']
    if mat['E'] is None:
        raise ValueError("No *Elastic (E) found")
    material_json = {"type": mat_type, "E": mat['E'], "rho": mat.get('rho', 0.0)}
    if mat_type == 'bar':
        material_json["area"] = data['thk'] if data['thk'] else 1.0
    else:
        material_json["nu"] = mat['nu'] if mat['nu'] is not None else 0.0
        if dim == 2:
            material_json["thk"] = data['thk'] if data['thk'] else 1.0

    # --- elements: KEEP ORIGINAL ABAQUS IDS (no renumber) ---
    elem_ids = sorted(data['elements'].keys())
    elements_json = []
    for eid in elem_ids:
        conn = data['elements'][eid]
        if len(conn) != nnode:
            raise ValueError("Element %d has %d nodes, expected %d for %s"
                             % (eid, len(conn), nnode, aba))
        elements_json.append({"id": eid, "connectivity": conn})

    group_json = {"type": elem_type, "material": 1, "elements": elements_json}

    # --- boundary conditions ---
    fixed = []
    prescribed = []

    def resolve_nodes(target):
        t = target.lower()
        if t in data['nsets']:
            return sorted(set(data['nsets'][t]))
        try:
            return [int(target)]
        except ValueError:
            print("  [warn] boundary/load target '%s' not found" % target)
            return []

    for row in data['boundary']:
        if not row:
            continue
        node_ids_bc = resolve_nodes(row[0])
        if not node_ids_bc:
            continue
        second = row[1] if len(row) > 1 else ''
        up = second.upper()
        if up in BC_TYPES:
            dofs = [d for d in BC_TYPES[up] if d <= 3]
            fixed.append({"nodes": node_ids_bc, "dof": [DOF_NAME[d] for d in dofs]})
        else:
            try:
                d1 = int(row[1])
                d2 = int(row[2]) if len(row) > 2 else d1
                val = float(row[3]) if len(row) > 3 else 0.0
            except (ValueError, IndexError):
                print("  [warn] unrecognized boundary row: %s" % row)
                continue
            dofs = [d for d in range(d1, d2 + 1) if d <= 3]
            if val == 0.0:
                fixed.append({"nodes": node_ids_bc,
                              "dof": [DOF_NAME[d] for d in dofs]})
            else:
                for n in node_ids_bc:
                    for d in dofs:
                        prescribed.append({"node": n, "dof": DOF_NAME[d], "value": val})

    # --- concentrated loads ---
    concentrated = []
    for target, dof, value in data['cload']:
        for n in resolve_nodes(target):
            concentrated.append({"node": n, "dof": DOF_NAME[dof], "value": value})

    # --- surface loads (TRVEC only) ---
    surface = []
    for surf_name, ltype, mag, direction in data['dsload']:
        if ltype != 'TRVEC':
            print("  [warn] surface load type '%s' on '%s' not supported "
                  "(only TRVEC); skipped. Use surface traction in ABAQUS."
                  % (ltype, surf_name))
            continue
        # normalize direction, scale by magnitude
        norm = math.sqrt(sum(d * d for d in direction))
        if norm < 1e-14:
            print("  [warn] zero direction vector on surface %s" % surf_name)
            continue
        comps = [mag * d / norm for d in direction]   # (fx, fy, fz)
        faces = resolve_surface_faces(data, surf_name)
        if not faces:
            print("  [warn] surface '%s' resolves to no element faces" % surf_name)
        for (eid, face_1based) in faces:
            for axis in range(dim):        # only x/y (2D) or x/y/z (3D)
                v = comps[axis]
                if abs(v) < 1e-14:
                    continue
                surface.append({
                    "element": eid,
                    "face": face_1based,
                    "dof": DOF_NAME[axis + 1],
                    "value": v,
                })

    model = {
        "title": "converted from %s" % aba,
        "dimension": dim,
        "nodes": nodes_json,
        "materials": [material_json],
        "element_groups": [group_json],
        "boundary_conditions": {"fixed": fixed, "prescribed": prescribed},
        "loads": {
            "concentrated": concentrated,
            "surface": surface,
            "body_force": data['body_force'] if data['body_force']
            else ([0.0, 0.0, 0.0] if dim == 3 else [0.0, 0.0]),
        },
    }
    return model


# --------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python inp2json.py <in.inp> [out.json]")
        sys.exit(1)
    inp = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else inp.rsplit('.', 1)[0] + '.json'

    print("Reading %s" % inp)
    data = parse_inp(inp)
    print("  element type: %s" % data['elem_type'])
    print("  nodes: %d, elements: %d" % (len(data['nodes']), len(data['elements'])))
    print("  nsets: %s" % sorted(data['nsets'].keys()))
    print("  elsets: %s" % sorted(data['elsets'].keys()))
    print("  surfaces: %s" % sorted(data['surfaces'].keys()))
    print("  boundary rows: %d, cload rows: %d, dsload rows: %d"
          % (len(data['boundary']), len(data['cload']), len(data['dsload'])))

    model = build_json(data)

    with open(out, 'w') as f:
        json.dump(model, f, indent=2)
    print("Wrote %s" % out)
    print("  fixed: %d, prescribed: %d, concentrated: %d, surface: %d"
          % (len(model['boundary_conditions']['fixed']),
             len(model['boundary_conditions']['prescribed']),
             len(model['loads']['concentrated']),
             len(model['loads']['surface'])))


if __name__ == '__main__':
    main()