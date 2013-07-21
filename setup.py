from distutils.core import setup
import sys

sys.path.append('kmc')
import kmc

sys.path.append('recoloring')
import recoloring

setup(name = 'recoloring',
      version = '0.1.0',
      author = 'Kris Gryte',
      author_email = 'kgryte@gmail.com',
      url = 'http://kgryte.com',
      download_url = '',
      description = 'Recoloring simulation.',
      long_description = recoloring.Recoloring.__doc__,
      package_dir = {'': 'recoloring-simulation'},
      py_modules = ['recoloring', 'kmc'],
      provides = ['recoloring', 'kmc'],
      keywords = 'kmc recoloring',
      license = 'MIT',
      classifiers = ['Natural Language :: English',
                     'Operating System :: OS Independent',
                     'Programming Language :: Python']
    )