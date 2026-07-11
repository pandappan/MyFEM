#!/usr/bin/env python3
"""
Convert ABAQUS .inp to STAP++ .dat format.
Specifically tailored for the plate-with-hole case.
"""
import re


def parse_inp(path):
    """Parse an ABAQUS .inp file. Returns nodes, elements, nsets, esets."""
    nodes = {}     # id -> (x, y)
    elements = {}  # id -> (n1, n2, n3, n4)
    nsets = {}     # name -> set of node ids
    esets = {}     # name -> set of element ids

    with open(path, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # ---- *Node ----
        if line == '*Node':
            i += 1
            while i < len(lines) and not lines[i].startswith('*'):
                parts = [p.strip() for p in lines[i].split(',')]
                if len(parts) >= 3:
                    nodes[int(parts[0])] = (float(parts[1]), float(parts[2]))
                i += 1
            continue

        # ---- *Element ----
        if line.startswith('*Element'):
            i += 1
            while i < len(lines) and not lines[i].startswith('*'):
                parts = [p.strip() for p in lines[i].split(',')]
                if len(parts) >= 5:
                    elements[int(parts[0])] = tuple(int(p) for p in parts[1:5])
                i += 1
            continue

        # ---- *Nset ----
        if line.startswith('*Nset'):
            m = re.search(r'nset=([^,]+)', line)
            name = m.group(1) if m else 'unknown'
            gen = 'generate' in line.lower()
            ids = []
            i += 1
            while i < len(lines) and not lines[i].startswith('*'):
                parts = [p.strip() for p in lines[i].split(',') if p.strip()]
                if gen:
                    s, e, st = int(parts[0]), int(parts[1]), int(parts[2])
                    ids += list(range(s, e + 1, st))
                else:
                    ids += [int(p) for p in parts]
                i += 1
            nsets[name] = set(ids)
            continue

        # ---- *Elset ----
        if line.startswith('*Elset'):
            m = re.search(r'elset=([^,]+)', line)
            name = m.group(1) if m else 'unknown'
            gen = 'generate' in line.lower()
            ids = []
            i += 1
            while i < len(lines) and not lines[i].startswith('*'):
                parts = [p.strip() for p in lines[i].split(',') if p.strip()]
                if gen:
                    s, e, st = int(parts[0]), int(parts[1]), int(parts[2])
                    ids += list(range(s, e + 1, st))
                else:
                    ids += [int(p) for p in parts]
                i += 1
            esets[name] = set(ids)
            continue

        i += 1

    return nodes, elements, nsets, esets


def write_dat(path, nodes, elements,
              left_nodes, left_bottom_nodes, right_elems,
              *, E=207000.0, nu=0.312, thk=1.0, traction=100.0):

    with open(path, 'w') as f:
        # ---- Header ----
        f.write("Plate with hole, uniaxial tension\n")
        f.write(f"2 {len(nodes)} 1 1\n")   # 2D, numnp, numeg=1, modex=1

        # ---- Nodes ----
        for nid in sorted(nodes):
            x, y = nodes[nid]
            if nid in left_bottom_nodes:
                bx, by = 1, 1
            elif nid in left_nodes:
                bx, by = 1, 1
            else:
                bx, by = 0, 0
            f.write(f"{nid} {bx} {by} {x:.9f} {y:.9f}\n")

        # ---- Element group header: etype=3(Q4_PS), NUME, NUMMAT ----
        f.write(f"3 {len(elements)} 1\n")

        # ---- Material set 1: nset rho E nu thk ----
        f.write(f"1 0.0 {E} {nu} {thk}\n")

        # ---- Elements: id n1 n2 n3 n4 mset ----
        for eid in sorted(elements):
            n = elements[eid]
            f.write(f"{eid} {n[0]} {n[1]} {n[2]} {n[3]} 1\n")

        # ---- Concentrated loads (empty) ----
        f.write("*CLOAD\n0\n")

        # ---- Prescribed displacements (empty) ----
        f.write("*PREDISPLACEMENT\n0\n")

        # ---- Surface loads on right edge ----
        # ABAQUS S4 face -> our 1-based face 4
        # dof=0 (UX), value=traction
        f.write("*SLOAD\n")
        f.write(f"{len(right_elems)}\n")
        for eid in sorted(right_elems):
            f.write(f"{eid} 4 1 {traction}\n")


if __name__ == '__main__':
    INP  = 'abaqus.inp'
    DAT  = 'MyFEM.dat'

    nodes, elements, nsets, esets = parse_inp(INP)
    left       = nsets.get('Part-1-1_left', set())
    left_bot   = nsets.get('Part-1-1_left_bottom', set())
    right_elems = esets.get('_Part-1-1_right_S4', set())

    print(f"Parsed: {len(nodes)} nodes, {len(elements)} elements")
    print(f"Left BC nodes: {len(left)}, Left-bottom pin: {len(left_bot)}")
    print(f"Right surface elements: {len(right_elems)}")

    write_dat(DAT, nodes, elements, left, left_bot, right_elems)
    print(f"Wrote {DAT}")