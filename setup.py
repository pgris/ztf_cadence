from setuptools import setup

# get the version here
pkg_vars = {}

with open("version.py") as fp:
    exec(fp.read(), pkg_vars)

setup(
    name='ztf_cadence',
    version=pkg_vars['__version__'],
    description='a set of scripts to perform cadence studies',
    url='http://github.com/pgris/ztf_cadence',
    author='Ph.Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['ztf_cadence', 'ztf_cadence_input',
              'ztf_metrics', 'ztf_cadence_utils'],
    package_data={'ztf_cadence_input': ['*.txt']},
    python_requires='>=3.5',
    zip_safe=False,
    install_requires=[
        'ztf_pipeutils>=0.1',
        'healpy>=1.13.0'
    ],
)
