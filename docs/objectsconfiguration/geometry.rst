Geometry
--------
The MOM6 geometry is defined within **MOM_input**, but a generic sea-ice geometry needs to be defined in a geometry section of a yaml file following the example below:

.. code-block:: yaml

     num_ice_cat: 5  # Number of sea ice categories
     num_ice_lev: 4  # Number of vertical levels within each seaice bins
     num_sno_lev: 1  #  "                      "             snow     "
