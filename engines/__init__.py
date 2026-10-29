"""Separate finite molecular and periodic engines behind one dispatcher."""


def calculate(config, directory, progress=print):
    if config.system == 'molecule':
        if config.calculation == 'optimize':
            from .optimizer import optimize_molecule
            return optimize_molecule(config, directory, progress)
        if config.calculation == 'frequencies':
            from .frequencies import frequency_molecule
            return frequency_molecule(config, directory, progress)
        from .pyscf_engine import calculate_molecule
        return calculate_molecule(config, directory, progress)
    from .periodic import calculate_periodic
    return calculate_periodic(config, directory, progress)
