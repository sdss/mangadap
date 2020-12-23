
Development scripts
-------------------

Directory holds a series of development scripts.

 - `build_bias_database.py`: Pull a bunch of data from MAPS files.

 - `select_representative_spectra.py`: Select representative MaNGA
   spectra from the full list based on some MAPS properties.

 - `generate_test_spectra.py`: Using the results of
   `select_representative_spectra.py`, compile the selected spectra into
   a single fits file for quick access.

 - `plot_selected.py`: Using the results of `generate_test_spectra.py`
   and `select_representative_spectra.py`, plot the spectra and models
   in a PDF document with one spectrum per page.

 -


