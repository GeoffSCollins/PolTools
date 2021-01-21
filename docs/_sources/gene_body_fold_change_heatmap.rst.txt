##############################
*Gene Body Fold Change Heatmap*
##############################
The ``gene_body_fold_change_heatmap`` tool generates a red/blue heatmap plotting the log 2 fold change of position and
quantity of 3' end sorted by gene length reads from two numerator datasets and two denominator datasets.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage
===============================
**Usage**:
::

  GC_bioinfo gene_body_fold_change_heatmap [-h] [-u upstream_distance]
                                                [-d distance_past_tes]
                                                [-b bp_width] [-w width]
                                                [-e height] [-m max_log2_fc]
                                                [--minor_ticks minor_ticks]
                                                [--major_ticks major_ticks]
                                                [-t [threads]]
                                                truQuant_output_file
                                                numerator_correction_factor_one
                                                numerator_seq_file_one
                                                numerator_correction_factor_two
                                                numerator_seq_file_two
                                                denominator_correction_factor_one
                                                denominator_seq_file_one
                                                denominator_correction_factor_two
                                                denominator_seq_file_two
                                                output_prefix

===========================================    =========================================================================================================================================================
Required Arguments                             Description
===========================================    =========================================================================================================================================================
**truQuant Output File**                       File ending in -truQuant_output.txt generated from `truQuant <https://geoffscollins.github.io/GC_bioinfo/truQuant.html>`_
**Numerator Correction Factor One**            Correction factor applied to the seq file data.
**Numerator Sequencing File One**              Bed formatted file from a sequencing experiment.
**Numerator Correction Factor One**            Correction factor applied to the seq file data.
**Numerator Sequencing File One**              Bed formatted file from a sequencing experiment.
**Denominator Correction Factor One**          Correction factor applied to the seq file data.
**Denominator Sequencing File One**            Bed formatted file from a sequencing experiment.
**Denominator Correction Factor One**          Correction factor applied to the seq file data.
**Denominator Sequencing File One**            Bed formatted file from a sequencing experiment.
**Output Prefix**                              Output filename will begin with the output prefix and also contain the run parameters and ends in gene_body_heatmap.tiff.
===========================================    =========================================================================================================================================================


===========================    ===============================================================================================================================================================
Optional Arguments             Description
===========================    ===============================================================================================================================================================
**-u, --upstream_distance**    Distance upstream of the max TSS from `truQuant <https://geoffscollins.github.io/GC_bioinfo/truQuant.html>`_ to show on the heatmap. Default is 50,000 bp.
**-d, --distance_past_tes**    Distance past the TES to show on the heatmap. Default is 50,000 bp.
**-b, --bp_width**             Total distance shown on the heatmap. Default is 400,000 bp.
**-w, --width**                Width of the heatmap in pixels. Default is 2,000 px.
**-e, --height**               Height of the heatmap in pixels. Default is 2,000 px.
**-g, --gamma**                Gamma correction of the heatmap. Default is 2.2, which is no gamma correction.
**-m, --max_log2_fc**            Maximum value to consider as black. Default is the max value found. Decreasing this number will make the image darker.
**--minor_ticks**              Distance in bp between minor tick marks. Default is 10,000 bp.
**--major_ticks**              Distance in bp between major tick marks. Default is 50,000 bp.
**-t, --threads**              Maximum number of threads. Default is the number of threads on the system. This program will not use more four threads.
===========================    ===============================================================================================================================================================


==========================================================================
Behavior
==========================================================================
``gene_body_fold_change_heatmap`` will create a file starting with the output prefix which contains the run parameters and ends in gene_body_heatmap.tiff.

For example:

.. image:: images/treatment_max_2.0_width_400000bp_gene_body_fold_change_heatmap.png

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

  $ GC_bioinfo gene_body_fold_change_heatmap seq_file-truQuant_output.txt 1 flavo.bed 1 flavo.bed 1 dmso.bed 1 dmso.bed treatment -m 2