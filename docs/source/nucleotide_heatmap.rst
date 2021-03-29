##############################
*Nucleotide Heatmap*
##############################
The ``nucleotide_heatmap`` tool creates heatmaps plotting the presence of nucleotides for each gene centered on the max TSS.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage
===============================
**Usage**:
::

  PolTools nucleotide_heatmap [-h]
                                max_tss_file region_width heatmap_width
                                vertical_average


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Max TSS FIle**               Bed formatted file containing the max TSSs you want to use to make the heatmaps. This file can be generated from the
                               `make_regions_file_centered_on_max_tss program <https://geoffscollins.github.io/PolTools/make_regions_file_centered_on_max_tss.html>`_ with a region size parameter of 1.
**Region width**               Number of base pairs to show on the heatmaps (will go upstream and downstream by width / 2).
**Heatmap Width**              Number of pixels wide the heatmaps will be.
**Vertical Average**           Without vertical averaging, the heatmaps will be tall and skinny. To make the image shorter, increasing the vertical average. This will average the pixels
                               vertically while making the height of the image (number of max TSSs / vertical average) tall.
===========================    =========================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``nucleotide_heatmap`` will create four files (one for each nucleotide -- A, T, G, C) named as the max TSS file with the
vertical averaging parameter and the nucleotide being plotted.

For example:

.. image:: images/tQ_max_tss_average_2_A.png

.. code-block:: bash

  $ head tQ_max_tss.bed
  chr1    959255  959256  NOC2L   46      -
  chr1    960632  960633  KLHL17  27      +
  chr1    966521  966522  PLEKHN1 8       +
  chr1    1000096 1000097 HES4    87      -
  chr1    1000295 1000296 ISG15   12      +
  chr1    1020119 1020120 AGRN    35      +
  chr1    1074306 1074307 RNF223  10      -
  chr1    1116106 1116107 C1orf159        9       -
  chr1    1231971 1231972 SDF4    321     -
  chr1    1232242 1232243 B3GALT6 174     +

  $ PolTools nucleotide_heatmap tQ_max_tss.bed 100 2000 2