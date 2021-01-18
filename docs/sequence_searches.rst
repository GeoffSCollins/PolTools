##############################
*Sequence Searches*
##############################
The ``sequence_searches`` tool is used determine the locations where a certain sequence (or multiple sequences) are found.


.. note::

  This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed and hg38.fa to be downloaded to
  ~/GC_bioinfo/static. This can be done by running the commands at the bottom of this page.

===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo sequence_searches <regions file> <Sequence,startPosition:endPosition>


======================================   =========================================================================================================================================================
Option                                   Description
======================================   =========================================================================================================================================================
**Regions Filename**                     Bed formatted file containing all the genes to quantify (regions will be determined
                                         from the 3' end of each region in this file.
**Sequence,startPosition:endPosition**   The sequence to search for (includes non standard bases like W) in the region
                                         (endPosition is inclusive).
======================================   =========================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``sequence_searches`` will output the regions file provided with new columns for the searching sequences.

For example:

.. code-block:: bash

  $ cat POLR2A_inr.bed
  chr17   7484355 7484375 POLR2A  0       +

  $ cat CCNT1_inr.bed
  chr12   48716696        48716716        CCNT1   0       -

  $ GC_bioinfo sequence_searches POLR2A_inr.bed GCTGCT,-3:3
  Chromosome      Left    Right   Gene    Score   Strand  GCTGCT
  chr17   7484355 7484375 POLR2A  0       +       True

===============================
Download hg38.fa
===============================
**Download hg38.fa**:
To download the hg38.fa (fasta file for the whole genome), run the following commands in the static directory:

.. code-block:: bash

  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz

