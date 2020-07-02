##############################
*Read Through Transcription*
##############################
The ``read_through_transcription`` tool computes the coverage of 3' ends of sequencing data around the 3' end of
features provided. This can be used to look at read through or runaway transcription.


.. note::

    This tool requires bedtools to be installed.

===============================
Usage and option summary
===============================
**Usage**:
::

  python3 read_through_transcription <Regions Filename> <TSR Filename> <Output Filename> \
                                          <Upstream Distance> <Downstream Distance> <Sequencing Files>


===========================    =========================================================================================================================================================
Option                         Description
===========================    =========================================================================================================================================================
**Regions Filename**           Bed formatted file containing all the genes to quantify (regions will be determined from the 3' end of each region in this file.
<<<<<<< HEAD
**TSR Filename**               `tsrFinder <https://github.com/P-TEFb/tsrFinderM1>`_ output file which will be blacklisted.
=======
**TSR Filename**               tsrFinder (https://github.com/P-TEFb/tsrFinderM1) output file which will be blacklisted.
>>>>>>> ff57e0a8b656a594f41d720c5b13546aa6a020be
                               Simply type *no* to not blacklist TSRs.
**Output Filename**            Name of the output file.
**Upstream Distance**          The number of base pairs to subtract from the left position.
**Downstream Distance**        The number of base pairs to add from the left position.
**Sequencing Files**           Sequencing files to quantify separated by spaces.
===========================    =========================================================================================================================================================



==========================================================================
Behavior
==========================================================================
``read_through_transcription`` will report the position relative to the 3' end of the regions provided and the sum
of the 3' reads at that position (50 bp intervals).

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

  $ python3 read_through_transcription.py gene_body_regions.bed no output.txt 100 500 control.bed treated.bed
  $ cat output.txt
  Position        control.bed treated.bed
  -100    19627   14509
  -50     17838   14471
  0       20637   16395
  50      24370   18634
  100     24902   18578
  150     25444   18615
  200     26958   19160
  250     28065   19877
  300     29358   20841
  350     27988   19534
  400     27068   19176
  450     26912   19786

