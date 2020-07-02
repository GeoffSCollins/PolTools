##############################
*truQuant*
##############################
The ``truQuant`` tool is used to build an annotation of the transcribed genes from PRO-Seq data. It also quantifies the
pause regions and gene bodies from the generated annotation while blacklisting enhancers, downstream promoters, and
non RNA Polymerase II transcripts.

.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ and `tsrFinder <https://github.com/P-TEFb/tsrFinderM1>`_ be installed.

===============================
Usage and option summary
===============================
**Usage**:
::

  python3 truQuant <Sequencing Files>


===========================    =========================================================================================================================================================
Option                         Description
===========================    =========================================================================================================================================================
**Sequencing Files**           Sequencing files to quantify separated by spaces.
===========================    =========================================================================================================================================================

.. note::

  The first sequencing file provided will be used to generate the annotation.

==========================================================================
Behavior
==========================================================================
``truQuant`` will generate search regions 1000 bp upstream of the 5' end of protein coding genes from GENCODE v32. Then,
tsrFinder will be run to determine the max TSR in the search region. Inside this TSR, the max TSS will be chosen as the
annotated 5' end. The pause region will be the 150 bp region surrounding the weighted average TSS (avgTSS) and gene body
will be the end of the pause region to the TES. TSRs in the gene with more than 30% of the reads as the max TSR will be
blacklisted. 5' ends in the pause regions will be quantified and 3' ends in the gene bodies will be quantified.

.. image:: images/truQuant.png

For example:

.. code-block:: bash

  $ head -n 5 control.bed
  chr1    10080   10380   K00294:149:H35VNBBXY:6:2108:3742:16524  255     -
  chr1    10563   10611   K00294:149:H35VNBBXY:6:1126:31730:23241 255     -
  chr1    10563   10600   K00294:149:H35VNBBXY:6:2206:29630:38627 255     -
  chr1    10564   10620   K00294:149:H35VNBBXY:6:1212:19441:27971 255     -
  chr1    10564   10611   K00294:149:H35VNBBXY:6:1211:31121:35022 255     -

  $ head -n 5 treated.bed
  chr1    10156   10374   K00294:149:H35VNBBXY:6:1209:23104:15891 255     -
  chr1    10564   10593   K00294:149:H35VNBBXY:6:2220:4888:19777  255     -
  chr1    10564   10597   K00294:149:H35VNBBXY:6:2104:25570:41651 255     -
  chr1    10564   10600   K00294:149:H35VNBBXY:6:1205:16407:42724 255     -
  chr1    10565   10597   K00294:149:H35VNBBXY:6:2221:2077:20709  255     -

  $ python3 truQuant.py control.bed treated.bed
  $ head -n 1 control-truQuant_output.txt # Printing the headers
  Gene    Chromosome      Pause Region Left       Pause Region Right      Strand  Total 5' Reads  MaxTSS  MaxTSS 5' Reads
  Weighted Pause Region Center    STDEV of TSSs   Gene Body Left  Gene Body Right Gene Body Distance
  control.bed Pause Region   treated.bed Pause Region   combined.bed Gene Body   treated.bed Gene Body

