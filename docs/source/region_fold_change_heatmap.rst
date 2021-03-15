##############################
*Region Fold Change Heatmap*
##############################
The ``region_fold_change_heatmap`` tool generates a log2 fold change heatmap of [5'/3'/pileup] reads from a series of sequencing datasets around provided regions.


.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo region_fold_change_heatmap [-h] [-m max_log2_fc]
                                             [-r repeat_amount]
                                             [-v vertical_averaging]
                                             [-t [threads]]
                                             [--minor_ticks minor_ticks]
                                             [--major_ticks major_ticks]
                                             read type regions_file
                                             output_prefix
                                             numerator_filename_one
                                             numerator_norm_factor_one
                                             numerator_filename_two
                                             numerator_norm_factor_two
                                             denominator_filename_one
                                             denominator_norm_factor_one
                                             denominator_filename_two
                                             denominator_norm_factor_two


====================================    =========================================================================================================================================================
Required Arguments                      Description
====================================    =========================================================================================================================================================
**Read Type**                           Either five, three, or whole corresponding to 5', 3', or pileup reads.
**Regions File**                        Bed formatted file containing regions of the same width.
**Numerator Sequencing File One**       Bed formatted file from a sequencing experiment.
**Numerator Norm Factor One**           Correction factor applied to the seq file data.
**Numerator Sequencing File Two**       Bed formatted file from a sequencing experiment.
**Numerator Norm Factor Two**           Correction factor applied to the seq file data.
**Denominator Sequencing File One**     Bed formatted file from a sequencing experiment.
**Denominator Norm Factor One**         Correction factor applied to the seq file data.
**Denominator Sequencing File Two**     Bed formatted file from a sequencing experiment.
**Denominator Norm Factor Two**         Correction factor applied to the seq file data.
**Output Prefix**                       Output filename will begin with the output prefix and also contain the run parameters and ends in region_heatmap.tiff.
====================================    =========================================================================================================================================================


=============================    ===============================================================================================================================================================
Optional Arguments               Description
=============================    ===============================================================================================================================================================
**-m, --max_log2_fc**            Maximum value to consider as black. Default is the max value found. Decreasing this number will make the image darker.
**-r, --repeat_amount**          Number of pixels that should be displayed for each base pair
**-v, --vertical_averaging**     Number of lines to average vertically.
**-g, --gamma**                  Gamma correction of the heatmap. Default is 2.2, which is no gamma correction.
**--minor_ticks**                Distance in bp between minor tick marks. Default is no ticks.
**--major_ticks**                Distance in bp between major tick marks. Default is no ticks.
=============================    ===============================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``region_fold_change_heatmap`` will generate a heatmap of

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

  $ head numerator_seq_file_one.bed
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

  $ GC_bioinfo region_heatmap five regions_centered_on_max_tss.bed numerator_seq_file_one.bed 1 numerator_seq_file_two.bed 1.15 \
    denominator_seq_file_one.bed 0.95 denominator_seq_file_two.bed 1.02 five_heatmap -r 20 -m 10 --minor_ticks 10 --major_ticks 50
