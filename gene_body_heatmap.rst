##############################
*Gene Body Heatmap*
##############################
The ``gene_body_heatmap`` tool computes the most commonly used transcript length for each region in the region file centered on the max TSS.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage
===============================
**Usage**:
::

  GC_bioinfo gene_body_heatmap [-h] [-u upstream_distance]
                                    [-d distance_past_tes] [-b bp_width]
                                    [-w width] [-e height] [-g gamma]
                                    [-m max_black] [--minor_ticks minor_ticks]
                                    [--major_ticks major_ticks]
                                    truQuant_output_file correction_factor
                                    seq_file output_prefix


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**truQuant Output File**       File ending in -truQuant_output.txt generated from `truQUant <truQuant.rst>`_
**Correction Factor**          Correction factor applied to the seq file data.
**Sequencing File**            Bed formatted file from a sequencing experiment.
**Output Prefix**              Output filename .# TODO
===========================    =========================================================================================================================================================


===========================    ===============================================================================================================================================================
Optional Arguments             Description
===========================    ===============================================================================================================================================================
**-u, --upstream_distance**    Distance upstream of the max TSS from `truQUant <truQuant.rst>`_ to show on the heatmap. Default is 50,000 bp.
**-d, --distance_past_tes**    Distance past the TES to show on the heatmap. Default is 50,000 bp.
**-b, --bp_width**             Total distance shown on the heatmap. Default is 400,000 bp.
**-w, --width**                Width of the heatmap in pixels. Default is 2,000 bp.
**-e, --height**               Height of the heatmap in pixels
**-g, --gamma**                Gamma correction of the heatmap.
**-m, --max_black**
**--minor_ticks**
**--major_ticks**
===========================    ===============================================================================================================================================================


==========================================================================
Behavior
==========================================================================
``gene_body_heatmap`` will report the most commonly used transcript length from the maxTSS (regions file centered) for each gene

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

  $ GC_bioinfo gene_body_heatmap seq_file-truQuant_output.txt 1.00 seq_file.bed seq