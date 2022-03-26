##############################
*tsrPicker*
##############################
The ``tsrPicker`` tool identifies transcription start regions from PRO-Cap, PRO-Seq, and related sequencing experiments.

===============================
Usage
===============================
**Usage**:
::

  PolTools tsrPicker [-h] [-r [radius]] seq_file min_seq_depth


================================    =========================================================================================================================================================
Required Arguments                  Description
================================    =========================================================================================================================================================
**Sequencing File**                 Bed formatted file from a sequencing experiment.
**Min Seq Depth**                   The minimum number of 5' reads to be considered as a TSR.
================================    =========================================================================================================================================================


===========================    ===============================================================================================================================================================
Optional Arguments             Description
===========================    ===============================================================================================================================================================
**-r, --radius**               Number of base pairs to expand the TSR from. Default is 5. For example, a value of 5 generates 11 bp TSRs as 5 base pairs are added on each side.
===========================    ===============================================================================================================================================================


==========================================================================
Behavior
==========================================================================
``tsrPicker`` generates a file named the sequencing file plus the minimum sequencing depth and -TSR.bed.
This bed formatted file contains TSRs with the number of 5' ends at the max TSS. ``tsrPicker`` selects the
base with the largest number of 5' ends and generates a TSR of size 1 + (2 * radius) centered on that base. This process
is repeated with all bases that are not within a TSR that meet the minimum sequencing depth.

For example:

.. code-block:: bash

  $ head seq_file.bed
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

  $ PolTools tsrPicker seq_file.bed 200
  $ head seq_file_20_20_30_600-TSR.tab
  chr1    156216434       156216445       TSR740  88595   +
  chr1    149832651       149832662       TSR741  24039   +
  chr1    203305513       203305524       TSR742  18139   +
  chr1    146376801       146376812       TSR743  16896   +
  chr1    16740510        16740521        TSR744  16140   +
  chr1    149851056       149851067       TSR745  15735   +
  chr1    16895974        16895985        TSR746  15223   +
  chr1    28648594        28648605        TSR747  15157   +
  chr1    85580755        85580766        TSR748  14528   +
  chr1    148522595       148522606       TSR749  13905   +