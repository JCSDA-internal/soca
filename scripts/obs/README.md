# Download observations from public servers
The scripts in the current directory (soca/scripts/obs) can be used to download oceanographic observations from public servers.

## Getting Started
The scripts are written in bash. They could be executed to any machine with *nix OS which supports Bash shell. 

### Prerequisites
For the data hosted at the Physical Oceanography Distributed Active Archive Center (PO.DAAC), the access to the observations is a two-step process:

* The user has to register at https://urs.earthdata.nasa.gov. The password used for accessing the data can be retrieved from "Drive API Credentials" (https://podaac-tools.jpl.nasa.gov/drive/login)

* The credentials to be stored at the ~/.netrc of the local machine as follows:

```
touch ~/.netrc
chmod 600 ~/.netrc
echo 'machine podaac-tools.jpl.nasa.gov login USERNAME_from_Drive_API_Credentials password PASSWORD_from_Drive_API_Credentials' >> ~/.netrc
```

Similarly, for ocean color data hosted at the NASA's Ocean Biology Processing Group (OBPG/OB.DAAC), e.g. from MODIS-Aqua, data access requires users to login to the OceanColor Web's data access points using their Earthdata Login credentials:

* Users need to register at https://urs.earthdata.nasa.gov and record USERNAME and PASSWD as their Earthdata Login credentials.

* Users need to modify the following script before attempting to download MODIS-Aqua ocean color data:
  source.coastwatch_oc_viirs_modis.sh: replacing USERNAME and PASSWD with their Login credentials in the following line:
  wget -q -O - $source_dir |grep OC| wget --user=USERNAME --password=PASSWD --auth-no-challenge=on --base https://oceandata.sci.gsfc.nasa.gov/ -N --wait=0.5 --random-wait --force-html -i -
     
## How to use
The users have to modify the following three variables at the get_obs.sh, according to their choices:

```
date_start=YYYYMMDD
date_end=YYYYMMDD
output_path=/path/where/the/obs/will/be/saved/
```

## Additional Information

* When wget is called, the first call does not take into consideration the content of the ~/.netrc file. Therefore it prompts a 401 error, and it calls wget for a second time with the credentials from ~/.netrc and the observations are downloaded.

* If there are no data for a specific date, an error message will appear. This message can mislead the user, but the process does not fail. 

