##############################
*Base Distribution*
##############################
The ``base_distribution`` tool computes the average base composition at each position of the given region.


.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed and hg38.fa to be downloaded to
    ~/PolTools/static. This can be done by running the commands at the bottom of this page.


===============================
Usage
===============================
**Usage**:
::

  PolTools base_distribution [-h] regions_file


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Regions Filename**           Bed formatted file containing all the regions you want to average the sequences. This file can be obtained by using the
                               `make_regions_file_centered_on_max_tss program <https://geoffscollins.github.io/PolTools/make_regions_file_centered_on_max_tss.html>`_
===========================    =========================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``base_distribution`` will report the position relative to the middle of the region and the percent occupancy of each base.

.. image:: images/base_distribution.png

For example:

.. code-block:: bash

  $ head regions_centered_on_max_tss.bed
  chr1    959251  959261  NOC2L   46      -
  chr1    960627  960637  KLHL17  27      +
  chr1    966516  966526  PLEKHN1 8       +
  chr1    1000092 1000102 HES4    87      -
  chr1    1000290 1000300 ISG15   12      +
  chr1    1020114 1020124 AGRN    35      +
  chr1    1074302 1074312 RNF223  10      -
  chr1    1116102 1116112 C1orf159        9       -
  chr1    1231967 1231977 SDF4    321     -
  chr1    1232237 1232247 B3GALT6 174     +

  $ PolTools base_distribution regions_centered_on_max_tss.bed
  Position        A       T       G       C
  -5.0    0.14846637102734664     0.1808019216555802      0.3518107908351811      0.3189209164818921
  -4.0    0.15973762010347375     0.13470066518847007     0.2900960827790096      0.41546563192904656
  -3.0    0.1271249076127125      0.17756836659275685     0.432649667405765       0.2626570583887657
  -2.0    0.07843680709534367     0.23697339246119734     0.3852549889135255      0.29933481152993346
  -1.0    0.021433850702143386    0.2683850702143385      0.0725240206947524      0.6376570583887657
  1.0     0.5554323725055432      0.00803769401330377     0.35559866962305986     0.08093126385809313
  2.0     0.16038433111603842     0.2354028085735403      0.3685328898743533      0.235679970436068
  3.0     0.2156319290465632      0.2953621581670362      0.2557280118255728      0.2332779009608278
  4.0     0.13848854397634885     0.27023281596452325     0.33379526977087953     0.25748337028824836
  5.0     0.12601626016260162     0.23706577974870657     0.2983185513673319      0.33859940872135996


===============================
Download hg38.fa
===============================
**Download hg38.fa**:
To download the hg38.fa (fasta file for the whole genome), run the following commands in the static directory:

.. code-block:: bash

  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz

