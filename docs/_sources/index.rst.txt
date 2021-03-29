Welcome to PolTools!
======================================

PolTools is a collection of bioinformatic programs analyzing PRO-Seq, PRO-Cap, and similar next generation sequencing
experiments analyzing nascent RNAs.

It is primary written in python3, but utilizes R for plotting and C++ for computationally expensive analyses. These
tools were developed in the `Price Lab <https://price.lab.uiowa.edu>`_

Installation
############

Two options are provided for installation. Both assume you already have R and bedtools installed. If you do not have R installed, please
install it as described `here <https://cran.r-project.org/doc/manuals/r-release/R-admin.html>`_. If you do not have
bedtools installed, please install it as described `in the documentation <https://bedtools.readthedocs.io/en/latest/content/installation.html>`_

1) Install using pip and pypi
Inside a terminal, run the following commands to install the necessary files

.. code-block:: bash

  $ pip3 install PolTools
  $ source /etc/bash_completion.d/PolTools-completion.bash
  $ PolTools build

2) Install from the github source

.. code-block:: bash

  $ git clone https://github.com/GeoffSCollins/PolTools.git
  $ sudo python3 PolTools/setup.py install
  $ source /etc/bash_completion.d/PolTools-completion.bash
  $ PolTools build

.. toctree::
   :maxdepth: 1
   :caption: Table of Contents
   :glob:

   *