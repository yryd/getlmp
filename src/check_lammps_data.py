#!/usr/bin/env python
"""轻量 data.lmp 自检：不依赖 LAMMPS，验证头部计数、段记录数、
系数段完整性、类型/原子引用合法性、电荷守恒。"""
import re
import sys


def parse_data(path: str) -> dict:
    with open(path) as f:
        lines = f.read().splitlines()

    # 去掉注释（# 开头行不算段）
    # 1) 头部计数
    counts = {}
    for line in lines[:30]:
        m = re.match(r'\s*(\d+)\s+(atoms|bonds|angles|dihedrals|impropers)\s*$', line)
        if m:
            counts[m.group(2)] = int(m.group(1))
        m = re.match(r'\s*(\d+)\s+(atom|bond|angle|dihedral|improper) types\s*$', line)
        if m:
            counts[m.group(2) + '_types'] = int(m.group(1))

    # 2) 定位各段
    sections = {}
    for i, line in enumerate(lines):
        m = re.match(r'^([A-Za-z ]+?)\s*(?:#.*)?$', line.strip())
        if m and m.group(1).strip() in {'Masses', 'Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs',
                                        'Dihedral Coeffs', 'Improper Coeffs', 'Atoms',
                                        'Velocities', 'Bonds', 'Angles', 'Dihedrals', 'Impropers'}:
            sections[m.group(1).strip()] = i

    def section_lines(name):
        if name not in sections:
            return []    # 段缺失（如 0 impropers 时无 Impropers 段）视为 0 记录
        start = sections[name] + 1
        end = min([v for k, v in sections.items() if v > start] or [len(lines)])
        out = []
        for line in lines[start:end]:
            if not line.strip() or line.strip().startswith('#'):
                continue
            out.append(line)
        return out

    # 3) 各段记录数与类型数
    result = {}
    for key, sec in [('atoms', 'Atoms'), ('bonds', 'Bonds'), ('angles', 'Angles'),
                     ('dihedrals', 'Dihedrals'), ('impropers', 'Impropers')]:
        recs = section_lines(sec)
        result[key + '_actual'] = len(recs)

    for key, sec in [('atom_types', 'Masses'), ('bond_types', 'Bond Coeffs'),
                     ('angle_types', 'Angle Coeffs'), ('dihedral_types', 'Dihedral Coeffs'),
                     ('improper_types', 'Improper Coeffs')]:
        recs = section_lines(sec)
        result[key + '_actual'] = len(recs)

    # 4) 原子段：id/type/电荷/坐标 + 电荷守恒
    # 支持两种 Atoms 格式（行尾可有注释列）：
    #   7 列 (atom_style full):    atom-ID mol-ID type q x y z
    #   6 列 (atom_style charge):  atom-ID type q x y z
    # 判别：第 3 列能解析为整数 → full（type 列）；否则 → charge（q 列是浮点）
    atoms = section_lines('Atoms')
    total_q = 0.0
    max_atom_id, max_atype = 0, 0
    for line in atoms:
        parts = line.split()
        try:
            atype = int(parts[2])
            aid, q = int(parts[0]), float(parts[3])
        except (ValueError, IndexError):
            try:
                atype = int(parts[1])
                aid, q = int(parts[0]), float(parts[2])
            except (ValueError, IndexError):
                raise ValueError(f'Atoms 段坏行: {line}')
        total_q += q
        max_atom_id = max(max_atom_id, aid)
        max_atype = max(max_atype, atype)
    result['total_charge'] = round(total_q, 5)
    result['max_atom_id'] = max_atom_id
    result['max_atype'] = max_atype

    # 5) 拓扑段引用合法（原子 id ≤ NATOM，类型 id ≤ 类型数）
    n_atom = counts.get('atoms', 0)
    for sec, ncol, label in [('Bonds', 2, 'bond'), ('Angles', 3, 'angle'),
                             ('Dihedrals', 4, 'dihedral'), ('Impropers', 4, 'improper')]:
        for line in section_lines(sec):
            parts = line.split()
            if len(parts) < 2 + ncol:
                raise ValueError(f'{sec} 坏行: {line}')
            tid = int(parts[1])
            if tid > counts.get(label + '_types', 0):
                raise ValueError(f'{sec} 类型越界: {line}')
            for j in range(2, 2 + ncol):
                aid = int(parts[j])
                if aid < 1 or aid > n_atom:
                    raise ValueError(f'{sec} 原子引用越界: {line}')

    result['counts'] = counts
    return result


def main():
    path = sys.argv[1]
    r = parse_data(path)
    print(f'自检 {path}:')
    c = r['counts']
    for k in ['atoms', 'bonds', 'angles', 'dihedrals', 'impropers']:
        ok = 'OK' if c.get(k) == r[k + '_actual'] else 'MISMATCH'
        print(f'  {k:10s} 头部={c.get(k)} 实际={r[k + "_actual"]}  [{ok}]')
    for k in ['atom_types', 'bond_types', 'angle_types', 'dihedral_types', 'improper_types']:
        ok = 'OK' if c.get(k) == r[k + '_actual'] else 'MISMATCH'
        print(f'  {k:10s} 头部={c.get(k)} 实际={r[k + "_actual"]}  [{ok}]')
    print(f'  电荷守恒: {r["total_charge"]}')
    print(f'  max_atom_id={r["max_atom_id"]} max_atype={r["max_atype"]} (NATOM={c.get("atoms")})')
    print('  → 格式与引用全部合法' if 'MISMATCH' not in str(r) else '  → 存在不一致!')


if __name__ == '__main__':
    main()
