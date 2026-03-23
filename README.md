# smbpal: SMB for paleo ice sheet simulations

This is a simple module that takes insolation, temperature and precipitation 
as input, and returns the annual SMB of an ice-covered domain,
using either the PDD or ITM melt models. 

## Dependencies

NetCDF. You know what to do. E.g.

```bash
brew install netcdf-fortran
```

## Getting started

1. Download the code.
2. Configure the Makefile (basically point to where the NetCDF libraries can be found).
3. Compile the static library and test program.
4. Run the test program.

```bash
git clone git@github.com:alex-robinson/smbpal.git
cd smbpal
# Edit existing config file in config/ or use one available. Then...
python config.py config/macbook_gfortran
make clean
make lib # Make static library that can be linked to in a program
make test # Make a test program that links to libsmbpal
./libsmbpal/bin/test_smbpal.x # Run the test SMB program test/test_smbpal.f90
```

That's it. You now have a working SMB module that will let you use the PDD or ITM 1-layer snowpack models.

## Under the hood

What if I want to use pieces of the library to build my own interface. Great! All physics subroutines/functions use native Fortran types (arrays/scalars/etc.) and are called for an individual point (independent of a larger grid). So you can use your own methods to manage things like parallelization of the domain. Insolation is calculated internally as needed (for ITM, not PDD) and the shortwave radiation is available.

Some things to note:

- src/smb_itm.f90: individual subroutines/functions related to calculating a snowpack model based on the insolation-temperature melt (ITM) method.
- src/smb_pdd.f90: the same, but for the positive degree day (PDD) method.

- `smbpal_update_monthly`: Wrapper to calculate monthly and annual SMB and surface temperature from either PDD or ITM. Input fields are interpolated smoothly to daily values and the snowpack is calculated over the whole year. The routine expects gridded data currently.

- `smbpal_update_itm`: Routine to calculate the daily step of the snowpack evolution using the ITM model. 

- `calc_insol_day_pt(day,lat,time_bp,S0,day_year,fldr)`: Function to calculate the insolation for a given day of the year `day`, latitude `lat` and time before present `time_bp`. Also specify how many days are in the year, typically `day_year=360`, and the folder where the insolation input files are located, like `fldr=libs/insol/input`.
