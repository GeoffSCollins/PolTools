##############################
*TPS Distance Per Gene*
##############################
The ``tps_distance_per_gene`` tool computes the most commonly used transcript length for each region in the region file centered on the max TSS.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage
===============================
**Usage**:
::

  GC_bioinfo tps_distance_per_gene [-h] [-t [threads]]
                                        regions_filename sequencing_files
                                        [sequencing_files ...]


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Regions filename**           Bed formatted file containing all the regions you want to quantify (must be centered on +1 nt). This file can be generated from the
                               `make_regions_file_centered_on_max_tss program <make_regions_file_centered_on_max_tss.rst>`_
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
``tps_distance_per_gene`` will report the most commonly used transcript length from the maxTSS (regions file centered) for each gene

For example:

.. code-block:: bash

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

  $ GC_bioinfo tps_distance_per_gene regions_centered_on_max_tss.bed seq_file.bed > seq_file_tps_distance_per_gene.txt
  $ head seq_file_tps_distance_per_gene.txt
  Gene    seq_file.bed
  NOC2L   32
  KLHL17  42
  PLEKHN1 45
  HES4    16
  ISG15   44
  AGRN    33
  RNF223  35
  C1orf159        N/A
  SDF4    68