##############################
*Divergent Pileup Metaplot*
##############################
The ``divergent_pileup_metaplot`` tool computes the coverage of reads data around the center of features provided.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.


===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo divergent_pileup_metaplot [-h] [-t [threads]]
                                            regions_file sequencing_files
                                            [sequencing_files ...]


===========================    =========================================================================================================================================================
Required Arguments                         Description
===========================    =========================================================================================================================================================
**Regions File**               Bed formatted file containing all the regions to quantify (+1 nucleotide centered). These regions can be made from the `make_regions_file_centered_on_max_tss program <make_regions_file_centered_on_max_tss.rst>`_
**Sequencing Files**           Bed formatted file from a sequencing experiment.
===========================    =========================================================================================================================================================

===========================    ===============================================================================================================================================================
Optional Arguments             Description
===========================    ===============================================================================================================================================================
**-t, --threads**              Maximum number of threads. Default is the number of threads on the system. This program will not use more threads than the number of sequencing files provided
**Sequencing Files**           Additional sequencing files can be provided.
===========================    ===============================================================================================================================================================


==========================================================================
Behavior
==========================================================================
``divergent_pileup_metaplot`` will report the position relative to the center of the regions provided and the average
of the reads at that position.

For example:

.. image:: images/divergent_pileup_metaplot.png

\


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

  $ head seq_file.bed
  chr1    11981   12023   A00876:119:HW5F5DRXX:1:2168:2248:1407   255     -
  chr1    13099   13117   A00876:119:HW5F5DRXX:1:2203:31403:26757 255     -
  chr1    13356   13423   A00876:119:HW5F5DRXX:1:2151:15808:7827  255     -
  chr1    13435   13477   A00876:119:HW5F5DRXX:1:2273:15781:19241 255     -
  chr1    13739   13772   A00876:119:HW5F5DRXX:1:2256:29966:10520 255     -
  chr1    13741   13773   A00876:119:HW5F5DRXX:1:2235:4101:11882  255     -
  chr1    14178   14203   A00876:119:HW5F5DRXX:1:2115:8241:31422  255     -
  chr1    14734   14768   A00876:119:HW5F5DRXX:1:2165:23764:2440  255     -
  chr1    14988   15012   A00876:119:HW5F5DRXX:1:2219:16134:32784 255     -
  chr1    18337   18362   A00876:119:HW5F5DRXX:1:2149:32054:31328 255     -

  $ GC_bioinfo divergent_pileup_metaplot regions_centered_on_max_tss.bed seq_file.bed
  Position       seq_file.bed sense strand       seq_file.bed divergent strand
  -5.0    108.42036215816704      -0.6227827050997783
  -4.0    114.89458610495196      -0.6073540280857354
  -3.0    125.8349963045085       -0.5766814486326681
  -2.0    138.23946784922396      -0.5629157427937915
  -1.0    147.7911123429416       -0.5778824833702882
  1.0     268.01856984478934      -0.5244826311899483
  2.0     277.60643015521066      -0.53640059127864
  3.0     289.13239098300073      -0.5285476718403548
  4.0     300.0270694752402       -0.5216186252771619
  5.0     306.0120103473762       -0.5210643015521065
