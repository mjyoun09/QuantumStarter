"""Verified basis/pseudopotential capabilities, readable without PySCF on Windows."""
import json
from pathlib import Path

CATALOG = json.loads((Path(__file__).resolve().parent / 'data' / 'element_capabilities.json').read_text(encoding='utf-8'))
MOLECULAR_CAPABILITIES = CATALOG['molecular']
PERIODIC_CAPABILITIES = CATALOG['periodic']


def molecular_core_counts(basis):
    return MOLECULAR_CAPABILITIES[basis]['ecp_core_electrons']


def periodic_supported_elements(basis):
    return frozenset(PERIODIC_CAPABILITIES[basis]['valence_electrons'])
