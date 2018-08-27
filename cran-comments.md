## Test environments
* OS X 10.13.5 (on local installation), R 3.5.1
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0
* win_builder
* x86_64-w64-mingw32 (R dev version: 2018-08-05 r75062)
* CRAN auto-check servers 
    * r-devel-linux-x86_64-debian-gcc (R dev version: 2018-07-31 r75041)
    * r-devel-windows-ix86+x86_64 (R dev version: 2018-07-31 r75040)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 3 NOTEs:

(OS X)
* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Joshua Hewitt <joshua.hewitt@colostate.edu>’
  New submission

(Ubuntu)
* checking installed package size ... NOTE
  installed size is 18.9Mb
  sub-directories of 1Mb or more:
  data   3.5Mb
  libs  15.1Mb
  
  telefit uses the Eigen C++ template library, which inflates the installed package's size.  The
  data directory is also large because it includes model output used to run examples quickly.

(Win-devel, CRAN auto-check servers)
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Joshua Hewitt <joshua.hewitt@colostate.edu>'
  New submission
  Possibly mis-spelled words in DESCRIPTION:
  al (9:349)
  et (9:346)
  geostatistical (9:121)
  teleconnection (9:84)
  
  These words are not misspelled.
