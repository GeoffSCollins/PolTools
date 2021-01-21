##############################
*Sequence Searches*
##############################
The ``sequence_searches`` tool is used determine if a sequence is present in a certain location around regions.


.. note::

  This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed and hg38.fa to be downloaded to
  ~/GC_bioinfo/static. This can be done by running the commands at the bottom of this page.

===============================
Usage and option summary
===============================
**Usage**:
::

  GC_bioinfo GC_bioinfo sequence_searches [-h] regions_filename search [search ...]

===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Regions filename**           Bed formatted file containing all the regions you want to quantify (must be centered on +1 nt). This file can be generated from the
                               `make_regions_file_centered_on_max_tss program <https://geoffscollins.github.io/GC_bioinfo/make_regions_file_centered_on_max_tss.html>`_
**Search**                     Search region and sequence formatted as follows: (Sequence),(-/+)left:(-/+)right. Ex: TATA,-30:-20
===========================    =========================================================================================================================================================

===========================    =========================================================================================================================================================
Optional Arguments             Description
===========================    =========================================================================================================================================================
**Search**                     Additional searches can be provided.
===========================    =========================================================================================================================================================



==========================================================================
Behavior
==========================================================================
``sequence_searches`` will output the regions file provided with new columns for the searches containing True or False values.

For example:

.. code-block:: bash

  $ head POLR2A_inr.bed
  chr17   7484355 7484375 POLR2A  0       +

  $ GC_bioinfo sequence_searches POLR2A_inr.bed GCTGC,-3:3
  Chromosome      Left    Right   Gene    Score   Strand  GCTGC
  chr17   7484355 7484375 POLR2A  0       +       True

===============================
Download hg38.fa
===============================
**Download hg38.fa**:
To download the hg38.fa (fasta file for the whole genome), run the following commands in the static directory:

.. code-block:: bash

  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz

