from setuptools import setup
setup(name='xbpy',
        version='0.1',
        description='Utility library intended to encapsulate common tasks used in molecular analysis.',
        url='',
        author='Finn Mier',
        license='MIT',
        dependencies=['numpy', 'rdkit', 'scipy', 'seaborn', 'tqdm'],
        packages=['xbpy.rdutil', 'xbpy.dispatch', 'xbpy.mathutils', 'xbpy.morgan', 'xbpy.interactions', 'xbpy.pymolutil', 'xbpy'],
        zip_safe=False)
