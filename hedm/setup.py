try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

config = {
	'description': 'Tools for working with high energy ' \
		'diffraction microscopy (HEDM)',
	'author': 'Branden Kappes',
	'url': 'URL at which to get it.',
	'download_url': 'Where to download it.',
	'author_email': 'bkappes@mines.edu',
	'version': '0.1',
	'install_requires': [
		'nose',
		'matplotlib',
		'numpy'],
	'packages': [
		'hedm',
		'hedm.plot'],
	'scripts': [],
	'name': 'hedm',
}
setup(**config)
