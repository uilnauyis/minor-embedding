Getting started with the D-Wave(TM) Quantum Computer Solver API Python Package
==============================================================================

Version 3.0.1

Installation requires setuptools; it is available from:

http://pypi.python.org/pypi/setuptools

Install the Solver API package by running the install.py script:

python install.py

If the installation fails, you may have downloaded the package for a different
platform or Python version.  To verify that the installation worked correctly,
start Python and type:

>>> import dwave_sapi2

You can also check that you can connect to the D-Wave(TM) Quantum Computer
Solver API.  You will need two pieces of information: the SAPI URL and an
authentication token.  The SAPI URL is listed on the "Solver API" page of the
web user interface.  Authentication tokens are also obtained from the web
user interface: click on "API Tokens" in the menu under your user name.  Test
the connection using these commands:

>>> import dwave_sapi2.remote
>>> conn = dwave_sapi2.remote.RemoteConnection(url, token)

If an error occurs, it is likely one of parameters is incorrect.  Otherwise,
you are ready to start using the D-Wave(TM) Quantum Computer Solver API.

You can find some example code in the examples subdirectory.


Linux users: if you see an error message like this:

  libssl3.so: cannot open shared object file: No such file or directory

your system is likely missing NSS libraries.  Install the libnss3 (Ubuntu,
Debian) or nss (Fedora, Red Hat) package.
