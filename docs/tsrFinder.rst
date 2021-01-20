##############################
*tsrFinder*
##############################
The ``tsrFinder`` tool identifies transcription start regions from PRO-Cap, PRO-Seq, and related sequencing experiments.
This tsrFinder is a much quicker version (with few output modifications) of the previous two tsrFinders: `one <https://github.com/P-TEFb/tsrFinder>`_ and
`two <https://github.com/P-TEFb/tsrFinderM1>`_.

===============================
Usage
===============================
**Usage**:
::

  GC_bioinfo tsrFinder [-h] [-t [threads]]
                            seq_file window_size min_seq_depth
                            min_avg_transcript_length max_fragment_size
                            chrom_size_file


================================    =========================================================================================================================================================
Required Arguments                  Description
================================    =========================================================================================================================================================
**Sequencing File**                 Bed formatted file from a sequencing experiment.
**Window Size**                     The size of the TSRs to find.
**Min Seq Depth**                   The minimum number of 5' reads to be considered as a TSR.
**Min avg Transcript Length**       The minimum average transcript length will eliminate TSRs from sequencing artifacts.
**Max fragment size**               The maximum transcript length for a read to be included in tsrFinder analysis.
**Chrom size file**                 A file containing the chromosome sizes. This can be optained using `fetchChromSizes <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>`_
                                    The hg38 chromsome size file can be found in the GC_bioinfo static directory.
================================    =========================================================================================================================================================


===========================    ===============================================================================================================================================================
Optional Arguments             Description
===========================    ===============================================================================================================================================================
**-t, --threads**              Maximum number of threads. Default is the number of threads on the system. This program will not use more threads than twice the number of chromosomes in the
                               chrom size file.
===========================    ===============================================================================================================================================================


==========================================================================
Behavior
==========================================================================
``tsrFinder`` generate a file named the sequencing file plus the tsrFinder parameters and -TSR.bed. This file contains
the chromosome, left position, right position, the read sum, number of 5' ends, strand, TSS left and right positions,
the TSS strength, and the position of the average TSS (weighted mean location).

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

  $ head hg38.chrom.sizes
  chr1    248956422
  chr2    242193529
  chr3    198295559
  chr4    190214555
  chr5    181538259
  chr6    170805979
  chr7    159345973
  chrX    156040895
  chr8    145138636
  chr9    138394717

  $ GC_bioinfo tsrFinder seq_file.bed 20 20 30 600 hg38.chrom.sizes
  $ head seq_file_20_20_30_600-TSR.tab
  chr1    629421  629441  738     20      +       629431  629432  4       629431
  chr1    629490  629510  2263    64      +       629494  629495  14      629500
  chr1    629564  629584  13877   273     +       629571  629572  188     629572
  chr1    629685  629705  939     21      +       629698  629699  6       629697
  chr1    629708  629728  1701    43      +       629723  629724  12      629719
  chr1    629740  629760  1911    63      +       629759  629760  21      629753
  chr1    629919  629939  1219    32      +       629929  629930  6       629931
  chr1    630666  630686  2277    44      +       630681  630682  21      630679
  chr1    630824  630844  658     21      +       630828  630829  4       630834
  chr1    630879  630899  2237    49      +       630893  630894  16      630890