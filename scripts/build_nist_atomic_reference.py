"""Build data/nist_atomic_reference.json from the NIST Atomic Reference Data pages (LDA and LSD columns).

Source: S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles, C. W. Clark,
"Atomic Reference Data for Electronic Structure Calculations", NIST,
https://math.nist.gov/DFTdata/atomdata/tables/ptable.html (energies in hartree).

Usage (pages are treated as untrusted text and only parsed with fixed patterns):
  python -I scripts/build_nist_atomic_reference.py --download <empty dir>     # fetch Z=1..36, 1 s apart
  python -I scripts/build_nist_atomic_reference.py --html-dir <dir with 01H.html ...>

Conventions taken from the NIST notation page: for LSD the first eigenvalue of each (nl) is the
majority spin ("down" in NIST wording) and the second, printed on the following line under the LSD
column, is the minority spin. Only the neutral-atom block of each page is kept; RLDA/ScRLDA are not.
Every number is parsed from the page text; nothing is typed by hand.
"""
import argparse
from datetime import date
import hashlib
import html
import json
from pathlib import Path
import re
import time
import urllib.request

SYMBOLS = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
           'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
BASE_URL = 'https://math.nist.gov/DFTdata/atomdata/tables/'
ENERGY_KEYS = ('Etot', 'Ekin', 'Ecoul', 'Eenuc', 'Exc')
NUMBER = r'-?\d+\.\d+'


def page_name(z):
    return f'{z:02d}{SYMBOLS[z - 1]}.html'


def download(directory, max_z):
    directory.mkdir(parents=True, exist_ok=True)
    for z in range(1, max_z + 1):
        request = urllib.request.Request(BASE_URL + page_name(z),
                                         headers={'User-Agent': 'Mozilla/5.0 (compatible; research-script)'})
        with urllib.request.urlopen(request, timeout=30) as response:
            (directory / page_name(z)).write_bytes(response.read())
        time.sleep(1.0)


def neutral_block(text, z):
    match = re.search(r'<PRE>(.*?)</PRE>', text, flags=re.S | re.I)
    if not match:
        raise ValueError(f'Z={z}: <PRE> block not found')
    plain = html.unescape(re.sub(r'<[^>]+>', '', match.group(1)))
    lines = plain.splitlines()
    header = next(i for i, line in enumerate(lines) if line.strip())
    head = re.fullmatch(r'\s*(\d+)\s+([A-Z][a-z]?)\s+(.*?)\s*', lines[header])
    if not head or int(head.group(1)) != z or head.group(2) != SYMBOLS[z - 1]:
        raise ValueError(f'Z={z}: unexpected header {lines[header]!r}')
    end = next((i for i in range(header + 1, len(lines)) if re.match(r'\s*\d+\s+[A-Z][a-z]?\+', lines[i])),
               len(lines))
    return head.group(3).replace('^', ''), lines[header + 1:end]


def parse_block(lines, z):
    energies = {'LDA': {}, 'LSD': {}}
    eigen_lda, eigen_lsd = {}, {}
    index = 0
    while index < len(lines):
        line = lines[index]
        energy = re.fullmatch(rf'\s*(\w+)\s*=\s*({NUMBER})\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*', line)
        orbital = re.fullmatch(rf'\s*(\d[spdf])\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*', line)
        if energy and energy.group(1) in ENERGY_KEYS:
            energies['LDA'][energy.group(1)] = float(energy.group(2))
            energies['LSD'][energy.group(1)] = float(energy.group(3))
        elif orbital:
            label = orbital.group(1)
            eigen_lda[label] = float(orbital.group(2))
            lsd_end = orbital.end(3)
            following = lines[index + 1] if index + 1 < len(lines) else ''
            second = re.match(rf'\s*({NUMBER})', following)
            if not second or second.end(1) != lsd_end:
                raise ValueError(f'Z={z} {label}: LSD minority value not aligned under the LSD column')
            eigen_lsd[label] = {'majority': float(orbital.group(3)), 'minority': float(second.group(1))}
            index += 1
        index += 1
    for column in ('LDA', 'LSD'):
        if set(energies[column]) != set(ENERGY_KEYS):
            raise ValueError(f'Z={z}: missing energy terms in {column}')
    if not eigen_lda or set(eigen_lda) != set(eigen_lsd):
        raise ValueError(f'Z={z}: orbital table incomplete')
    return {'LDA': {**energies['LDA'], 'eigenvalues': eigen_lda},
            'LSD': {**energies['LSD'], 'eigenvalues': eigen_lsd}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--html-dir', type=Path)
    source.add_argument('--download', type=Path, help='empty directory to download the pages into first')
    parser.add_argument('--max-z', type=int, default=36)
    parser.add_argument('--accessed', default=date.today().isoformat(), help='date the pages were fetched')
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'data' / 'nist_atomic_reference.json')
    args = parser.parse_args()
    directory = args.download or args.html_dir
    if args.download:
        download(directory, args.max_z)
    elements = {}
    for z in range(1, args.max_z + 1):
        raw = (directory / page_name(z)).read_bytes()
        configuration, lines = neutral_block(raw.decode('utf-8', errors='replace'), z)
        elements[SYMBOLS[z - 1]] = {'Z': z, 'configuration': configuration, **parse_block(lines, z),
                                    'page': BASE_URL + page_name(z), 'page_sha256': hashlib.sha256(raw).hexdigest()}
    data = {
        'source': {'title': 'Atomic Reference Data for Electronic Structure Calculations',
                   'authors': 'S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles, C. W. Clark',
                   'publisher': 'National Institute of Standards and Technology (NIST)',
                   'index': BASE_URL + 'ptable.html', 'accessed': args.accessed, 'units': 'hartree',
                   'columns': {'LDA': 'nonrelativistic spin-unpolarized LDA (Vosko-Wilk-Nusair correlation)',
                               'LSD': 'nonrelativistic local spin density'},
                   'lsd_eigenvalue_order': 'majority spin first (NIST "down"), minority spin second',
                   'energy_terms': {'Etot': 'total', 'Ekin': 'kinetic', 'Ecoul': 'Hartree electron-electron',
                                    'Eenuc': 'electron-nucleus', 'Exc': 'exchange-correlation'},
                   'parser': 'scripts/build_nist_atomic_reference.py (numbers parsed from the page text)'},
        'elements': elements}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes((json.dumps(data, indent=1) + '\n').encode('utf-8'))
    print(f'{len(elements)} elements -> {args.output}')


if __name__ == '__main__':
    main()
