from sys import exit, version_info
from os.path import dirname, realpath, join
from setuptools.command import easy_install
from pkg_resources import get_build_platform

if version_info[:2] != (2, 7):
    exit("Python 2.7 required")

plat = get_build_platform().rsplit('-', 1)[-1]

eggdir = join(dirname(realpath(__file__)), 'eggs')
eggs = (
    "sapiremote-3.0.1-py2.7-linux-{plat}.egg",
    "sapilocal-3.0.1-py2.7-linux-{plat}.egg",
    "fix_variables-3.0.1-py2.7-linux-{plat}.egg",
    "find_embedding-3.0.1-py2.7-linux-{plat}.egg",
    "qsage-3.0.1-py2.7-linux-{plat}.egg",
    "dwave_sapi2-3.0.1-py2.7.egg")

for egg in eggs:
    easy_install.main(['-H', '', join(eggdir, egg.format(plat=plat))])
