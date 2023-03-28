.. CryoDRGN documentation master file, created by
   sphinx-quickstart on Mon Sep 12 11:12:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CryoDRGN's documentation!
====================================

Welcome to ❄️🐉  **cryoDRGN**'s documentation! This page contains user guides for the ``cryodrgn`` open source software package, including in-depth tutorials for training and analyzing cryoDRGN models.

Quick start
-----------

CryoDRGN can be installed with ``pip``.
We recommend installing ``cryodrgn`` in a separate anaconda environment.
See :doc:`/pages/installation` for more details and advanced installation instructions. ::

   $ pip install cryodrgn

The entry point for all cryoDRGN commands are the ``cryodrgn`` and ``cryodrgn_utils`` executables::

   (cryodrgn) $ cryodrgn -h
   (cryodrgn) $ cryodrgn_utils -h

Use the ``-h`` flag for help, e.g. to display all available subcommands or ``cryodrgn <command> -h`` to see the available parameters to each command.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pages/intro
   pages/installation
   pages/empiar_tutorial
   pages/landscape_analysis
   pages/large_datasets
   pages/cryodrgn2
   pages/faq


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
