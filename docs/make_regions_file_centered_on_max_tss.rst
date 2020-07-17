#######################################
*Make Regions File Centered on Max TSS*
#######################################
The ``make_regions_file_centered_on_max_tss`` tool generates a bed formatted region file centered on the max TSS from truQuant


===============================
Usage and option summary
===============================
**Usage**:
::

  python3 make_regions_file_centered_on_max_tss.py  <truQuant output file> <region size>


===========================    =========================================================================================================================================================
Option                         Description
===========================    =========================================================================================================================================================
**truQuant output file**       `truQuant <https://github.com/GeoffSCollins/GC_bioinfo/blob/master/docs/truQuant.rst>`_ output file (not the paused, gene body, or blacklisted regions
                               file!
**Region Size**                Integer value of the number of base pairs the region will include
===========================    =========================================================================================================================================================

==========================================================================
Behavior
==========================================================================
``make_regions_file_centered_on_max_tss`` will report the generated regions centered on the max TSS from truQuant in bed format.
For example:

.. code-block:: bash

  $ grep "POLR2A" truQuant_output.bed
  POLR2A  chr17   7484299 7484449 +       5278    7484369 2478    7484372 12.970928230422533      7484450 7514618 30168   5277    1169    3667.5363333645746 812.4597259244244

  $ python3 make_regions_file_centered_on_max_tss.py CCNT1_inr.bed 20 | grep "POLR2A"
  chr17   7484359 7484379 POLR2A  2478    +