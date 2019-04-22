========================
NDEx TCGA Content Loader
========================


.. image:: https://img.shields.io/pypi/v/ndextcgaloader.svg
        :target: https://pypi.python.org/pypi/ndextcgaloader

.. image:: https://img.shields.io/travis/coleslaw481/ndextcgaloader.svg
        :target: https://travis-ci.org/coleslaw481/ndextcgaloader

.. image:: https://readthedocs.org/projects/ndextcgaloader/badge/?version=latest
        :target: https://ndextcgaloader.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Loads TCGA data into NDEx


* Free software: BSD license
* Documentation: https://ndextcgaloader.readthedocs.io.


Tools
-----

* **ndexloadtcga.py** -- Loads TCGA into NDEx_

Dependencies
------------

* `ndex2 3.1.0a1 <https://pypi.org/project/ndex2/3.1.0a1/>`_
* `ndexutil 0.2.0a1 <https://pypi.org/project/ndexutil/0.2.0a1/>`_

Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block::

   git clone https://github.com/coleslaw481/ndextcgaloader
   cd ndextcgaloader
   make dist
   pip install dist/ndextcgaloader*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image 

.. code-block::

   make

Output:

.. code-block::

   clean                remove all build, test, coverage and Python artifacts
   clean-build          remove build artifacts
   clean-pyc            remove Python file artifacts
   clean-test           remove test and coverage artifacts
   lint                 check style with flake8
   test                 run tests quickly with the default Python
   test-all             run tests on every Python version with tox
   coverage             check code coverage quickly with the default Python
   docs                 generate Sphinx HTML documentation, including API docs
   servedocs            compile the docs watching for changes
   testrelease          package and upload a TEST release
   release              package and upload a release
   dist                 builds source and wheel package
   install              install the package to the active Python's site-packages
   dockerbuild          build docker image and store in local repository
   dockerpush           push image to dockerhub


Configuration
-------------

The **ndexloadtcga.py** requires a configuration file in the following format be created.
The default path for this configuration is :code:`~/.ndexutils.conf` but can be overridden with
:code:`--conf` flag.

**Format of configuration file**

.. code-block::

    [<value in --profile (default ndextcgaloader)>]

    user = <NDEx username>
    password = <NDEx password>
    server = <NDEx server(omit http) ie public.ndexbio.org>

The NDEx UUID needed for **style** can be obtained by uploading the :code:`style.cx` file found under
the :code:`data/` directory of this repository. NOTE: The network needs to be uploaded to the same
server as defined in **style** :code:`public.ndexbio.org` is NDEx_ production. Also the network needs
to be visible to the **user**

**Example configuration file**

.. code-block::

    [ndextcgaloader_dev]
    user = joe123
    password = somepassword123
    server = dev.ndexbio.org

    [ndextcgaloader_prod]
    user = joe_p
    password = joes_unbreakable_password
    server = ndexbio.org

Needed files
------------

Three files are needed to run this script: **loadplan.json**, **networks.txt** and **style.cx**. These files are loacted in **data** directory.  They are specified with **--loadplan**, **--networklistfile** and **--style** command-line arguments, accordingly.
For example:
.. code-block::
   ndexloadtcga.py --loadplan ./data/loadplan.json --networklistfile ./data/networks.txt --style ./data/style.cx 

Usage
-----

For information invoke :code:`ndexloadtcga.py -h`

In addition to the three required files listed in the previous section, we need to configure and specify profile (production or development), and working directory where tsv networks will be created before uploading to the server.

The entire command is thus
.. code-block::
 ndexloadtcga.py --loadplan <loadplan> --networklistfile <networks file> --style <style> --profile <profile> -datadir <datadir>

**Example usage**

**TODO:** Add information about example usage

.. code-block::

   ndexloadtcga.py --loadplan ./data/loadplan.json --networklistfile ./data/networks.txt --style ./data/style.cx --profile ndextcgaloader_dev -datadir ./networks


   ndexloadtcga.py --loadplan ./data/loadplan.json --networklistfile ./data/networks.txt --style ./data/style.cx --profile ndextcgaloader_prod -datadir ./networks


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

**TODO:** Add information about example usage


.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` coleslawndex/ndextcgaloader:0.1.0 ndexloadtcga.py --conf conf # TODO Add other needed arguments here


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
