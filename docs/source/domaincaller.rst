Detect Single-level TAD
***********************
Domaincaller is a wrap of the first 2 steps of `HiTAD <https://xiaotaowang.github.io/TADLib/hitad.html>`_,
which can be used as a method to detect single-level domains from genome-wide contact matrix. Similar
to HiTAD, it also takes a `cool <https://github.com/mirnylab/cooler>`_ file as input, and then calculates
the adaptive DI track and performs Hidden Markov Model (HMM) to predict TAD boundary locations.


Usage
^^^^^
To run *domaincaller*, just follow the pseudo command below::

    $ domaincaller --uri /path/to/the/cool/URI -O test.tad.bed --DI-output test.DIs.bedGraph --removeCache
  
Type ``domaincaller`` or ``domaincaller -h`` on your terminal to print detailed help information for each parameter.






