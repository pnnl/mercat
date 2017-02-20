from setuptools import setup, find_packages

required = open('dependencies.txt').read().splitlines()
required = [l.strip() for l in required
            if l.strip() and not l.strip().startswith('#')]

setup(
    name='mercat',
    version="0.1",
    author='Ajay Panyala',
    author_email='ajay@example.com',
    url='https://github.com/pnnl/mercat',
    packages=['mercat'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'mercat=mercat.mercat:mercat_main',
        ],
    }
)


