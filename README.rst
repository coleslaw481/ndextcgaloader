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

Python application that loads TCGA networks into NDEx_

**1\)** downloads network files in text format from server specified by ``--dataurl`` argument (default is `https://github.com/iVis-at-Bilkent/pathway-mapper/tree/master/samples <https://github.com/iVis-at-Bilkent/pathway-mapper/tree/master/samples>`_)

**2\)** the list of files to be downloaded is specified by ``--networklistfille`` argument (default is ``networks.txt`` that comes with the distribution of this utility)

**3\)** the files are downloaded to a directory specified by ``--datadir`` argument (default is ``network`` in ndextcgaloader installation directory)

**4\)** after that, utility generates CX networks and uploads them to the server; below is a brief description of how downloaded text files is transformed to CX format:
 * text file is opened for reading and description of network is extracted (if it is there)
 * two panda dataframes created, one is initialized with the definitions of nodes and another one with the definitions of edges, read from the text file 
 * a dictionary, with node Ids as keys and node names as values, is initialized
 * all node names are then checked if they make valid HGNC names (if not, they are recorded to ``reports/invalid_protein_names.tsv``),  and if they are not nested (i.e., complex nodes containing another complex nodes, for example, FAMILY containing FAMILY; all nested nodes are recorded to ``reports/nested_nodes.tsv``)
 * nested nodes (if any) are normalized, i.e., if node type proteinfamily A has as its' member node type proteinfamily B, and node B has three genes (C, D, E), then node type proteinfamily B is removed, and genes C, D, E are made members of  proteinfamily A
 * duplicate edges, if any, are removed (leaving one edge), some of edge and node headers are renamed for readability, node and edge dataframes are joined into one dataframe
 * orphan gene nodes are removed 
 * members for complex nodes (node types other than genes) are generated
 * then the pandas dataframe is saved to tsv file
 * a network in NiceCX is generated from the panda dataframe (network descripiton extracted earlier is used); this network is saved on the disk
 * after saving, the network in CX is used to replace the existing network on the server, or upload to server (if network doesn't exist there)
    
**5\)** to connect to NDEx server and upload generated in CX format networks, a configuration file must be passed with ``--conf`` parameter. If ``--conf`` is not specified, the configuration ``~/{confname}`` is examined.

Dependencies
------------

* `ndex2 3.1.0a1 <https://pypi.org/project/ndex2/3.1.0a1/>`_
* `ndexutil 0.2.0a1 <https://pypi.org/project/ndexutil/0.2.0a1/>`_

Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block:: python

   git clone https://github.com/ndexcontent/ndextcgaloader.git
   cd ndextcgaloader
   make dist
   pip install dist/ndextcgaloader*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image

.. code-block:: python

   make

Output:

.. code-block:: python

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

The ``ndexloadtcga.py`` requires a configuration file in the following format be created.
The default path for this configuration is ``~/.ndexutils.conf`` but can be overridden with
``--conf`` flag.

**Format of configuration file**

.. code-block:: python

    [<value in --profile (default ndextcgaloader)>]
    user = <NDEx username>
    password = <NDEx password>
    server = <NDEx server(omit http), i.e., public.ndexbio.org>


**Example of a default configuration for Development server in ~/.ndexutils.conf:**

Default configuration is defined in the section ``[ndextcgaloader]`` :

.. code-block:: python

    [ndextcgaloader]
    user = joe123
    password = somepassword123
    server = dev.ndexbio.org


**Example of configuration for Production server in ~/.ndexutils.conf:**

.. code-block:: python

    [ndextcgaloader_prod]
    user = joe_p
    password = joes_unbreakable_password
    server = ndexbio.org

Usage
-----

**Running with default configuration**

To run utility with the above default config, it is suffice to call utility with no arguments:

.. code-block:: python

    ndexloadtcga.py

This will upload networks to account ``joe123`` on server ``dev.ndexbio.org`` (specified in ``[ndextcgaloader]`` section of ``~/.ndexutils.conf``)


**Running with explicitly specified configuration**

To make ``ndexloadtcga.py`` upload networks to account ``joe_p`` on ``ndexbio.org``:

.. code-block:: python

    ndexloadtcga.py --profile ndextcgaloader_prod


Needed files
------------

Three files needed to run this script are:

.. code-block:: python

   loadplan.json
   networks.txt
   style.cx

These files are located in NDEX TCGA Loader installation directory.  They are used by the script by default. Users, however, may want to specify their own loadplan, list of networks or style instead of the provided default ones. To do so, please use ``--loadplan``, ``--networklistfile`` and/or ``--style`` command-line arguments. For example, in order to use your own style defined in ``my_style.cx``:

.. code-block:: python

   ndexloadtcga.py --style my_style.cx


``reports`` directory
---------------------

``ndexloadtcga.py`` creates ``reports`` directory with two files in ``tsv`` format:

.. code-block:: python

   nested_nodes.tsv
   invalid_protein_names.tsv

``nested_nodes.tsv`` contains list of complex nodes (nodes that are not proteins) that have other complex nodes as members. ``invalid_protein_names.tsv`` contains list of invalid names found in networks.  These files are provided for information/debugging purpose and can be safely deleted.


More information
----------------

For more information invoke 

.. code-block:: python

   ndexloadtcga.py -h



Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

**TODO:** Add information about example usage


.. code-block:: python

   docker run -v `pwd`:`pwd` -w `pwd` coleslawndex/ndextcgaloader:0.1.0 ndexloadtcga.py --conf conf # TODO Add other needed arguments here


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
