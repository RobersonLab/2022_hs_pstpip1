from setuptools import setup

setup(
	name = "burdens_code",
	packages = ['burdens_code'],
	version = "v1.0.0",
	author = "David J. Morales-Heil",
	install_requires = ['numpy', ],
	entry_points = {'console_scripts':["burdens_code = burdens_code.__main__:run"]}
)