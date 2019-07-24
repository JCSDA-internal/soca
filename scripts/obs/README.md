# Download observations from public servers
The scripts in the current directory (soca/scripts/obs) can be used to download observations from public servers.

## Getting Started
The scripts are written in bash. They could be executed to any machine with linux OS.

### Prerequisites
For the data hosted at the Physical Oceanography Distributed Active Archive Center (PO.DAAC), the access to the observations is two step process:

* The user has to register at https://urs.earthdata.nasa.gov. The same credentials are used for the PO.DAAC.

* The credentials to be stored at the ~/.netrc of the local machine

```
touch ~/.netrc
chmod 600 ~/.netrc
echo 'machine podaac-tools.jpl.nasa.gov login USERNAME_from_PO.DAAC password PASSWORD_from_PO.DAAC' >> ~/.netrc

```

## How to use
The user has to modify the following three variables at the get_obs.sh, according to their choices


```
date_start= 
date_end=
output_path=

```
## Additional Information

* When wget is called, the first call does not take into consideration the content of the ~/.netrc file, therefore it prompts a 401 error and it calls wget for a second time with the credentials from ~/.netrc and the observations are downloaded.

* If there are no data for specific date, a error message will appear. This message can mislead the user. 

