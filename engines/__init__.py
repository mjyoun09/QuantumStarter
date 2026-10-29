"""Separate finite molecular and periodic engines behind one dispatcher."""


def calculate(config, directory, progress=print):
    if config.system == 'molecule':
        from .pyscf_engine import calculate_molecule
        return calculate_molecule(config, directory, progress)
    from .periodic import calculate_periodic
    return calculate_periodic(config, directory, progress)
