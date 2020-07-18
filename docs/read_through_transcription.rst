##############################
*Read Through Transcription*
##############################
The ``read_through_transcription`` tool computes the coverage of 3' ends of sequencing data around the 3' end of
features provided. This can be used to look at read through or runaway transcription.


.. note::

    This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo read_through_transcription.py <Regions Filename> <TSR Filename> \
          <Upstream Distance> <Downstream Distance> <Interval Distance> <Sequencing Files>


===========================    =========================================================================================================================================================
Option                         Description
===========================    =========================================================================================================================================================
**Regions Filename**           Bed formatted file containing all the genes to quantify (regions will be determined from the 3' end of each region in this file.
**TSR Filename**               `tsrFinder <https://github.com/P-TEFb/tsrFinderM1>`_ output file which will be blacklisted.
                               Simply type *no* to not blacklist TSRs.
**Upstream Distance**          The number of base pairs to subtract from the left position.
**Downstream Distance**        The number of base pairs to add from the left position.
**Interval Distance**          The size of sub-regions to split the regions into.
**Sequencing Files**           Sequencing files to quantify separated by spaces.
===========================    =========================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``read_through_transcription`` will report the position relative to the 3' end of the regions provided and the sum
of the 3' reads at that interval.

For example:


.. image:: images/read_through_transcription.png

\

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

  $ GC_bioinfo read_through_transcription.py gene_body_regions.bed no output.txt 100 500 50 control.bed treated.bed
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

