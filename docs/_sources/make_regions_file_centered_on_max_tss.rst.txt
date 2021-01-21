#######################################
*Make Regions File Centered on Max TSS*
#######################################
The ``make_regions_file_centered_on_max_tss`` tool generates a bed formatted region file centered on the max TSS from truQuant


===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo make_regions_file_centered_on_max_tss [-h]
                                                   truQuant_file
                                                   region_size


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**truQuant output file**       `truQuant <https://geoffscollins.github.io/GC_bioinfo/truQuant.html>`_ output file (not the paused, gene body, or blacklisted regions
                               file!)
**Region Size**                Integer value of the number of base pairs the region will include
===========================    =========================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``make_regions_file_centered_on_max_tss`` will report the generated regions centered on the max TSS from truQuant in bed format.
For example:

.. image:: images/make_regions_file_centered_on_max_tss.png

\

.. code-block:: bash

  $ head seq_file-truQuant_output.txt
  Gene    Chromosome      Pause Region Left       Pause Region Right      Strand  Total 5' Reads  MaxTSS  MaxTSS 5' Reads Weighted Pause Region Center    STDEV of TSSs   Gene Body Left  Gene Body Right Gene Body Distance      seq_file.bed Pause Region   seq_file.bed Gene Body
  NOC2L   chr1    959177  959327  -       194     959255  46      959250  13.306459171023036      944203  959177  14974   194     18
  KLHL17  chr1    960552  960702  +       234     960632  27      960626  25.417791063821863      960702  965719  5017    234     17
  PLEKHN1 chr1    966439  966589  +       25      966521  8       966513  19.47408534437497       966589  975865  9276    25      11
  HES4    chr1    1000013 1000163 -       239     1000096 87      1000086 27.14758979723915       998962  1000013 1051    239     68
  ISG15   chr1    1000204 1000354 +       160     1000295 12      1000278 36.24344768368484       1000354 1014540 14186   160     111
  AGRN    chr1    1020042 1020192 +       112     1020119 35      1020116 25.189637892253575      1020192 1056118 35926   112     76
  RNF223  chr1    1074208 1074358 -       32      1074306 10      1074284 32.567238138964136      1070967 1074208 3241    32      8
  C1orf159        chr1    1116028 1116178 -       51      1116106 9       1116103 19.81136532595448       1081818 1116028 34210   51      11
  SDF4    chr1    1231907 1232057 -       1105    1231971 321     1231978 23.701136922154493      1216908 1231907 14999   1097    177

  $ GC_bioinfo make_regions_file_centered_on_max_tss seq_file-truQuant_output.txt 1 > tQ_max_tss.bed

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