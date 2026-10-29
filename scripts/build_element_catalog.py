"""Regenerate the frontend capability table from the installed PySCF data."""
import json
from pathlib import Path
import sys
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular import ELEMENTS
import pyscf
from pyscf import gto
from pyscf.pbc.gto import basis as pbc_basis, pseudo

warnings.simplefilter('ignore', UserWarning)
names = ('sto-3g', '6-31g', '6-31g*', '6-31+g*', 'cc-pvdz',
         'def2-svp', 'def2-tzvp', 'dyall-v2z', 'dyall-v3z')
data = {'schema_version': 1, 'pyscf_version': pyscf.__version__, 'molecular': {}, 'periodic': {}}
for name in names:
    support, cores = [], {}
    for e in ELEMENTS:
        try:
            if not gto.basis.load(name, e):
                continue
        except (gto.basis.BasisNotFoundError, KeyError, AttributeError):
            continue
        support.append(e)
        if name.startswith('def2-'):
            ecp = gto.basis.load_ecp(name, e)
            if ecp:
                cores[e] = int(ecp[0])
    data['molecular'][name] = {'elements': support, 'ecp_core_electrons': cores}
for name in ('gth-dzvp-molopt-sr', 'gth-szv-molopt-sr', 'gth-dzvp', 'gth-tzvp'):
    valence = {}
    for e in ELEMENTS:
        try:
            basis = pbc_basis.load(name, e)
            potential = pseudo.load('gth-pbe', e)
            if basis and potential:
                valence[e] = int(sum(potential[0]))
        except (gto.basis.BasisNotFoundError, KeyError, AttributeError):
            pass
    data['periodic'][name] = {'pseudo': 'gth-pbe', 'valence_electrons': valence}
output = ROOT / 'data' / 'element_capabilities.json'
output.parent.mkdir(exist_ok=True)
output.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print({name: len(info['elements']) for name, info in data['molecular'].items()})
print({name: len(info['valence_electrons']) for name, info in data['periodic'].items()})
