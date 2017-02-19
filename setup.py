import io
from os.path import dirname, join
from setuptools import setup, find_packages


def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
    if "__version__" in line:
      if '"' in line:
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]

setup(
    name='mercat',
    version=get_version("mercat/__init__.py"),
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


