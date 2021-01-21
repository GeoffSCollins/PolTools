##############################
*Five Prime Metaplot*
##############################
The ``five_prime_metaplot`` tool computes the average coverage of 5' ends of sequencing data around the center of features provided.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo five_prime_metaplot [-h] [-t [threads]]
                                      regions_file sequencing_files
                                      [sequencing files ...]


===========================    =========================================================================================================================================================
Required Arguments                         Description
===========================    =========================================================================================================================================================
**Regions File**               Bed formatted file containing all the regions to quantify (+1 nucleotide centered). These regions can be made from the `make_regions_file_centered_on_max_tss program <make_regions_file_centered_on_max_tss.rst>`_
**Sequencing Files**           Bed formatted file from a sequencing experiment.
===========================    =========================================================================================================================================================


===========================    =========================================================================================================================================================
Optional Arguments                         Description
===========================    =========================================================================================================================================================
**-t, --threads**              Maximum number of threads. Default is the number of threads on the system. This program will not use more threads than the number of sequencing files provided
**Sequencing Files**           Additional sequencing files can be provided.
===========================    =========================================================================================================================================================


==========================================================================
Behavior
==========================================================================
``five_prime_metaplot`` will report the position relative to the center of the regions provided and the average
of the 5' reads at that position.

For example:

.. image:: images/five_prime_metaplot.png

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

  $ GC_bioinfo five_prime_metaplot regions_centered_on_max_tss.bed seq_file.bed
  Position        seq_file.bed 5' sense strand    seq_file.bed 5' divergent strand
  -5.0    6.177753141167775       -0.0291019955654102
  -4.0    7.360218033998522       -0.039634146341463415
  -3.0    11.54619364375462       -0.027439024390243903
  -2.0    13.114098300073909      -0.013396156688839615
  -1.0    10.78279748706578       -0.059959349593495935
  1.0     120.73272357723577      -0.0041574279379157425
  2.0     11.020140428677013      -0.020417590539541758
  3.0     12.267645971914265      -0.017738359201773836
  4.0     11.616962305986696      -0.01524390243902439
  5.0     7.43080192165558        -0.01681448632668145