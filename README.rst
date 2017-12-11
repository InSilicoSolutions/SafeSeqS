SafeSeqS – Safe Sequencing System
---------------------------------------------

Python package for processing of SafeSeqS data (https://www.ncbi.nlm.nih.gov/pubmed/21586637).   The pipeline analyzes multiple
samples in fastq sequencing files to identify mutations present in a small fraction of DNA templates.


Documentation
-------------

Please see the `documentation <https://github.com/InSilicoSolutions/SafeSeqS/wiki>`_ on github for complete usage details.

Installation
------------

Use pip to install the SafeSeqS package:

.. code-block:: bash

    $ pip install safeseqs

SafeSeqS depends on several other packages and these will be installed automatically by pip if they
are not already installed.  

**Required packages:**

* scipy
* pywin32 (for Windows machines)

Note pip is sometimes not able to install pywin32.  If so, please download the executable from https://sourceforge.net/projects/pywin32/files/pywin32/ that
matches your python version and operating system.  Then re-run pip install safeseqs.
