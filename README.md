# RADPROFfitter
Radial profile fitter for "A modular set of synthetic spectral energy distributions for YSOs" [Robitaille, T. P., 2017](https://ui.adsabs.harvard.edu/abs/2017A%26A...600A..11R/abstract)

This is python based module for observational data radial profile's fitting. Now it's support MSX, WISE, IRAS and custom image data of YSOs.

Before use it you need to check if your fits files contains correct header:
- correct WCS
- information about wavelength (for custom image it must be in microns)
- correct information to get pixel Scale
- Flux in Jy/pixel for custom image_data


If you discover bugs, have issues with the code or suggestions, please feel free to open an issue on Github.

## Acknowledgements

Radial profile fitter for a modular set of synthetic spectral energy distributions for YSOs were implemented with the financial support of grant
[18-72-10132](https://rscf.ru/en/project/18-72-10132/) of the [Russian Science
Foundation](https://rscf.ru/en/).


## References

- *HYPERION: An open-source parallelized three-dimensional dust continuum radiative transfer code*; [Robitaille, T. P., 2011](https://ui.adsabs.harvard.edu/abs/2011A%26A...536A..79R/abstract)

- *A modular set of synthetic spectral energy distributions for young stellar objects*; [Robitaille, T. P., 2017](https://ui.adsabs.harvard.edu/abs/2017A%26A...600A..11R/abstract)
