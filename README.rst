Module-painter: Fast assessment of modular recombinations in bacteriophages
===========================================================================

.. image:: https://img.shields.io/github/workflow/status/Puumanamana/module-painter/test?label=tests&logo=github
   :target: https://github.com/Puumanamana/module-painter/actions?query=workflow

Citation (Work in progress)
---------------------------
Cédric Arisdakessian, Mahdi Belcaid, Anne Bergeron and Guylaine Poisson
Module-painter: Fast assessment of modular recombinations in bacteriophages

Description
-----------
Module-painter is a tool for identifying the traces of a parent viral population on a target population. By identifying recombination patterns, it clusters the viruses into subpopulations based on their shared history.
Module-painter is composed of two main steps. First, we determine a minimal set of segments from the parent population covering each child’s genome. Second, we identify all recombinations shared by multiple children in order to cluster them into subpopulations.

Install
-------

Dependencies:

- `minimap2` (or `blastn` if `--aligner blastn` is set)
- python dependencies (see setup.py file)
- (optional) `pycairo` for graph plotting

.. code-block:: bash

   git clone https://github.com/Puumanamana/module-painter.git \
   && cd module-painter \
   && pip install .

Basic usage
-----------
To use module-painter on your own dataset, you need to provide two mandatory arguments:
- The viral population (list of fasta files)
- The regex patterns or keywords identifying the children population (unless the --rotate-parent flag is set, in which case the parent population rotates among all input fasta files (leave one out strategy)

Example:

.. code-block:: bash

   module-painter run -p data/*.fasta -c children

.. code-block:: bash

   module-painter run -p data/*.fasta --rotate-parent

For a comprehensive list of all options, see the command help:

.. code-block:: bash

   module-painter run -h

Contribute
----------

- Issue Tracker: `github <https://github.com/Puumanamana/module-painter/issues>`_
- Source Code: `github <https://github.com/Puumanamana/module-painter>`_
