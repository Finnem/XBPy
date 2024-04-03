from setuptools import setup
setup(name='xbpy',
        version='0.1',
        description='Utility library intended to encapsulate common tasks used in molecular analysis.',
        url='',
        author='Finn Mier',
        license='MIT',
        packages=['xbpy.rdutil', 'xbpy.dispatch', 'xbpy.math', 'xbpy.morgan', 'xbpy'],
        zip_safe=False)
