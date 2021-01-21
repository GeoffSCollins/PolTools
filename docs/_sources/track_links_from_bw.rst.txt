##############################
*Track Links From Bigwigs*
##############################
The ``track_links_from_bw`` tool prints the UCSC Genome browser links for the supplied bigwig files.

===============================
Usage
===============================
**Usage**:
::

  GC_bioinfo track_links_from_bw [-h] bigwig_files [bigwig_files ...]


===========================    =========================================================================================================================================================
Required Arguments             Description
===========================    =========================================================================================================================================================
**Bigwig files**               Bigwig formatted files
===========================    =========================================================================================================================================================



==========================================================================
Behavior
==========================================================================
``track_links_from_bw`` will print the UCSC genome browser links for the supplied bigwig files. Once these links are in a file,
open the file in a text editor to replace the <server> text with the server you are hosting the data.

For example:

.. code-block:: bash

  $ GC_bioinfo track_links_from_bw seq_file-FW.bw seq_file-RV.bw
  track type=bigWig visibility=full name='seq file FW' autoScale=on alwaysZero=on windowingFunction=maximum negateValues=off color=0,0,0 altColor=0,0,0 bigDataUrl=<server>/seq_file-FW.bw
  track type=bigWig visibility=full name='seq file RV' autoScale=on alwaysZero=on windowingFunction=maximum negateValues=on color=0,0,0 altColor=0,0,0 bigDataUrl=<server>/seq_file-RV.bw