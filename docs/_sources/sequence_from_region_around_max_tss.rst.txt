##############################
*Sequence from Around Max TSS*
##############################
The ``sequence_from_region_around_max_tss`` tool returns the sequence of selected regions around the maxTSS.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage
===============================
**Usage**:
::

  GC_bioinfo sequence_from_region_around_max_tss [-h] max_tss_file left right


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Max TSS FIle**               Bed formatted file containing the max TSSs you want to use to make the heatmaps. This file can be generated from the
                               `make_regions_file_centered_on_max_tss program <https://geoffscollins.github.io/GC_bioinfo/make_regions_file_centered_on_max_tss.html>`_ with a region size parameter of 1.
**Left**                       Left part of the region to get the sequence.
**Right**                      Right part of the region to get the sequence.
===========================    =========================================================================================================================================================


==========================================================================
Behavior
==========================================================================
``sequence_from_region_around_max_tss`` will print out a fasta formatted file containing the sequences.

For example:

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

  $ GC_bioinfo sequence_from_region_around_max_tss tQ_max_tss.bed -5 5 > output.tmp
  $ head output.tmp
  >NOC2L
  TGCACGCTTC
  >KLHL17
  CCGGCAGTCT
  >PLEKHN1
  TGTACGACTC
  >HES4
  CGCGCGCGGG
  >ISG15
  CCTCCGACAC
