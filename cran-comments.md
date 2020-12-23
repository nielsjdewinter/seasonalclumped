## New Submission
This is a new submission of the seasonalclumped package

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

There were 2 NOTES:

* Possibly mis-spelled words in DESCRIPTION:
    al (14:16)
    de (14:3)
    et (14:13)
    speleothems (13:55)

  These are parts of academic citations and cannot be replaced
  All remaining spelling mistakes flagged by "devtools::spell_check()" are jargon or names

* Found the following (possibly) invalid URLs:
  URL: http://doi.org/bvpzws (moved to https://doi.org/10.1016/0168-9622(86)90057-6)
    From: man/binning_seasonality.Rd
          man/carbmodel.Rd
          man/optimization_seasonality.Rd
          man/oxygen_isotope_seasonality.Rd
          man/smoothing_seasonality.Rd
    Status: 200
    Message: OK

  Checked this short DOI and it is valid and links to the correct reference.

## Downstream dependencies
There are currently no downstream dependencies for this package