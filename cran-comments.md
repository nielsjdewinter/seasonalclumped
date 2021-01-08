## New Submission
This is a resubmission. In this version I have made the following changes:
* Reduced package title length to below 45 characters
* Simplified examples to allow execution during check()
* Removed instances where local parameters were changed
* Added missing \arguments and \values to Rd files
* Corrected error in example of oxygen_isotope_seasonality()

## Test environments
* local Win 10 Pro install, R 4.0.2
* Windows Server 2008 R2 SP1 (R-devel on R-hub builder)
* Fedora Linux (R-devel on R-hub builder)
* Ubuntu 16.04 LTS (R-release on R-hub builder)
* Ubuntu 16.04 (on travis ci) R 4.0.2
* OS X 10.13 (on travis ci) R 4.0.3

## R CMD check results
There were no ERRORs, WARNINGSs or NOTEs

## R-hub Builder results
There were no ERRORs or WARNINGs.

There were 1 NOTE:

* Possibly mis-spelled words in DESCRIPTION:
    al (14:16)
    de (14:3)
    et (14:13)
    speleothems (13:55)

  These are parts of academic citations and cannot be replaced
  All remaining spelling mistakes flagged by "devtools::spell_check()" are jargon or names

## Downstream dependencies
There are currently no downstream dependencies for this package