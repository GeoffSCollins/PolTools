##############################
*Region Heatmap*
##############################
The ``region_heatmap`` tool generates a heatmap of [5'/3'/pileup] reads from a sequencing dataset around provided regions.


.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo region_heatmap [-h] [-m max_black] [-r repeat_amount]
                                 [-v vertical_averaging] [-g gamma]
                                 read type regions_file sequencing_file
                                 norm_factor output_prefix


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Read Type**                  Either five, three, or whole corresponding to 5', 3', or pileup reads.
**Regions File**               Bed formatted file containing regions of the same width.
**Sequencing File**            Bed formatted file from a sequencing experiment.
**Norm Factor**                Correction factor applied to the seq file data.
**Output Prefix**              Output filename will begin with the output prefix and also contain the run parameters and ends in region_heatmap.tiff.
===========================    =========================================================================================================================================================


=============================    ===============================================================================================================================================================
Optional Arguments               Description
=============================    ===============================================================================================================================================================
**-m, --max_black**              Maximum value to consider as black. Default is the max value found. Decreasing this number will make the image darker.
**-r, --repeat_amount**          Number of pixels that should be displayed for each base pair
**-v, --vertical_averaging**     Number of lines to average vertically.
**-g, --gamma**                  Gamma correction of the heatmap. Default is 2.2, which is no gamma correction.
**--minor_ticks**                Distance in bp between minor tick marks. Default is no ticks.
**--major_ticks**                Distance in bp between major tick marks. Default is no ticks.
=============================    ===============================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``region_heatmap`` will generate a heatmap of

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

  $ GC_bioinfo region_heatmap five regions_centered_on_max_tss.bed seq_file.bed 1 five_heatmap -r 20 -m 10\
    --minor_ticks 10 --major_ticks 50
