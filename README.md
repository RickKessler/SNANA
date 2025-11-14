# SNANA
Supernova Analysis package.

Public data/filters/models/etc are in SNDATA_ROOT (~2 GB):
   https://zenodo.org/records/12655677

Read documentation in SNANA/doc.

SNANA code moved from snana.uchicago.edu to Github (May 2019).

SNANA tutorial:
  https://kicp.uchicago.edu/~kessler/SNANA_Tutorial/SNANA_Tutorial_2023-05.pdf


If 'git pull' results in  "detected dubious ownership" error on a compute cluster, try the following:  
   git config --global --add safe.directory <SNANA_DIR>  
   git pull  


Using Autotools to build:
autoreconf -fi   # -f=force -i = install
./configure
make             # builds all
make install     # installs all *.exe in bin directory

to build one executable:
cd src
make snana.exe

to clean:
make clean     # removes *.o and *.exe files
make distclean # does make clean plus removes files created with ./configure


