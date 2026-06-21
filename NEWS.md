# dyadicMarkov 0.1.1

* Updated package wording and metadata for the CRAN submission.
* Added S3 classes and print/summary support for pattern and case identification results.
* Added S3 classes for empirical count matrices and MLE transition probability matrices while preserving
  ordinary matrix behavior.
* Added two synthetic 90-point example datasets for package workflow examples.
* Rewrote the workflow vignette around the built-in univariate and bivariate example datasets.
* Improved internal input validation for count, estimation, and pattern-identification functions.
* Updated tests and documentation for the new S3 return objects.
* Improved validation for extreme state-space inputs, non-finite chain values, and malformed empirical
  matrices.
* Refactored selected internal validation and AIC helper code to reduce function complexity while preserving exported
  behavior.
* Improved bivariate count validation coverage for unsupported and malformed inputs.

# dyadicMarkov 0.1.0

* Initial CRAN submission.
