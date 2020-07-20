##############################
*Initiator Reads*
##############################
The ``inr_reads`` tool computes the number of 5' end reads at the max TSS of a regions file.

===============================
Usage
===============================
**Usage**:
::

  GC_bioinfo inr_reads <regions file> <sequencing files>


===========================    =========================================================================================================================================================
Option                         Description
===========================    =========================================================================================================================================================
**Regions Filename**           Bed formatted file containing all the regions you want to quantify (must be centered on +1 nt)
**Sequencing Files**           Sequencing files to quantify separated by spaces.
===========================    =========================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``inr_reads`` will report the number of 5' ends at the +1 nt in the center of the regions file.

For example:

.. code-block:: bash

  $ head promoters.bed
  chr4    498989  499489  PIGG    38      +
  chr4    674291  674791  MYL5    17      +
  chr4    705581  706081  PCGF3   41      +
  chr4    932209  932709  TMEM175 33      +
  chr4    1011383 1011883 FGFRL1  29      +
  chr4    1289643 1290143 MAEA    83      +
  chr4    1346997 1347497 UVSSA   47      +
  chr4    1721268 1721768 TACC3   56      +
  chr4    1793042 1793542 FGFR3   216     +
  chr4    1871130 1871630 NSD2    11      +

  $ head control.bed
  chr1    11242   11265   A00876:65:HLHG7DRXX:2:2274:6922:19867   255     -
  chr1    11295   11325   A00876:65:HLHG7DRXX:1:2165:16432:10473  255     -
  chr1    11691   11829   A00876:65:HLHG7DRXX:1:2203:10059:23171  255     +
  chr1    11760   11914   A00876:65:HLHG7DRXX:1:2269:19560:5181   255     -
  chr1    12019   12114   A00876:65:HLHG7DRXX:2:2276:27208:3928   255     -
  chr1    12285   12312   A00876:65:HLHG7DRXX:2:2220:16767:29778  255     -
  chr1    12291   12376   A00876:65:HLHG7DRXX:1:2112:17400:17096  255     -
  chr1    12302   12342   A00876:65:HLHG7DRXX:2:2145:2483:31219   255     -
  chr1    12310   12348   A00876:65:HLHG7DRXX:1:2273:11478:4492   255     -
  chr1    12381   12401   A00876:65:HLHG7DRXX:1:2219:8829:14434   255     +

  $ GC_bioinfo inr_reads promoters.bed control.bed > output.tmp
  $ head output.tmp
  Gene    control.bed
  PIGG    38
  MYL5    17
  PCGF3   41
  TMEM175 33
  FGFRL1  29
  MAEA    83
  UVSSA   47
  TACC3   56
  FGFR3   216
