##############################
*Sequence Searches*
##############################
The ``sequence_searches`` tool is used determine if a sequence is present in a certain location around regions.


.. note::

  This tool requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

===============================
Usage and option summary
===============================
**Usage**:
::

  PolTools PolTools sequence_searches [-h] regions_filename search [search ...]

===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Regions filename**           Bed formatted file containing all the regions you want to quantify (must be centered on +1 nt). This file can be generated from the
                               `make_regions_file_centered_on_max_tss program <https://geoffscollins.github.io/PolTools/make_regions_file_centered_on_max_tss.html>`_
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

  $ PolTools sequence_searches POLR2A_inr.bed GCTGC,-3:3
  Chromosome      Left    Right   Gene    Score   Strand  GCTGC
  chr17   7484355 7484375 POLR2A  0       +       True
