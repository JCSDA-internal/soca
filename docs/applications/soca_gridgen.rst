.. _applications_soca_gridgen:

soca_gridgen.x
================

Description
--------------

   This application was created to bypass some of the MOM6 initialization when instantiating a :ref:`Geometry` object. It generates a grid using MOM6's grid generation tools and saves the parameters relevant to the Geometry. The following is saved in a fms compliant netcdf file, soca_gridspec.nc:

   - Longitude
   - Latitude
   - Cell area
   - Rossby radius
   - Ocean/Land Mask

Configuration
--------------

 The YAML configuration file has 2 sections:

.. code-block:: yaml

   geometry:

   gridgen:

The :ref:`geometry <Geometry>` section contains the few arguments needed to setup the sea-ice geometry.
The gridgen section is currently empty, and is a place holder for future configurations.
