#
#                               3DVAR (No FGAT)
#
#     <--------------------- window_length (ex:6 hours) --------------------->
# t_start                              t_mid                               t_end 
# ---|-----------|-----------|-----------|-----------|-----------|-----------|
#    |                                   |                                   |      
# window_begin: 2015-01-01T09:00:00Z     |                           2015-01-01T15:00:00Z
#                                        |
#                       bkg/ana/incr: 2015-01-01T12:00:00Z
#
#

import yaml
import datetime
from pydic2yml import *

def main():

    # Parse 3DVAR parameters
    with open('3dvar.yml','r') as ymlfilehandle:
        try:
            da_params = yaml.load(ymlfilehandle)
        except yaml.YAMLError as err:
            print err

    # Create yml file for soca/oops 3dvar
    #Start of DA window
    t_start=datetime.datetime.strptime(da_params.get('start_date'),'%Y-%m-%dT%H:%M:%SZ')

    print 'Start date:',t_start
    
    #DA window length [hours]
    dt=da_params.get('window_length')
    dt_mid=dt/2
    dt=datetime.timedelta(hours=dt)
    dt_mid=datetime.timedelta(hours=dt_mid)

    #Date at the middle of the DA window
    t_mid=t_start+dt_mid
    print 'Middle of window:',t_mid
    
    #Ocean file name (ex: MOM.res_Y2016_D003_S64800.nc)
    #day of year
    doy = t_mid.toordinal() - datetime.date(t_mid.year, 1, 1).toordinal() + 1
    #second of day
    sod = t_mid.hour*3600+t_mid.minute*60+t_mid.second
    ocnfname='MOM.res_Y'+str(t_mid.year)+'_D'+str(doy).zfill(3)+'_S'+str(sod).zfill(5)+'.nc'

    #Sea-ice file name (ex: 20160103.180000.ice_model.res.nc)
    icefname=str(t_mid.year)+str(t_mid.month).zfill(2)+str(t_mid.day).zfill(2)+ \
     '.'+str(t_mid.hour).zfill(2)+str(t_mid.minute).zfill(2)+ \
         str(t_mid.second).zfill(2)+'.ice_model.res.nc'

    model_vars=["cicen","hicen","hsnon","tsfcn","qsnon","sicnk","qicnk",
            "socn","tocn","ssh"]
    window_begin=t_start.isoformat()+'Z'
    window_length='PT'+str(dt.seconds/3600)+'H'


    #Number of outer loop
    nouter = da_params.get('number_of_outer_loops')


    #Setup dict for yml parsing
    model=Model()
    bkg_error=Covariance(date=t_mid.isoformat()+'Z')
    bkg_mid=BackgroundFile(date=t_mid.isoformat()+'Z',
                       basename='./RESTART/',
                       ocn_filename=ocnfname,
                       ice_filename=icefname)
    bkg_mid=[bkg_mid.todict()]

    #Observation types
    obstypes=ObsType(obstype='SeaIceFraction',
                     obsfileout='Data/test.out',
                     obsfilein='/scratch4/NCEPDEV/ocean/scrub/Guillaume.Vernieres/JEDI/seaice_obs.nc',
                     obsvalue='ObsVal',
                     covariance='diagonal',
                     obserror='ObsErr')
    obstypes=[obstypes.todict()]
    jo=Jo(obstypes)
    jb=Jb(backgrounds=bkg_mid,covariance=bkg_error)
    J=CostFunction(variables=model_vars,Jb=jb,Jo=jo,window_begin=window_begin,window_length=window_length)

    iterations=Iteration(tlm='IdTLM',resolution={},tstep=model.dt,ninner=1)

    tmpiter=[]
    for iter in range(nouter):
        tmpiter.append(iterations.todict())
    var=Variational(tmpiter)
        
    minimizer=Minimizer()

    output=Output()

    data = dict(
        test_framework_runtime_config = "--log_level=test_suite",
        resolution = dict(),
        model = model.todict(),
        cost_function = J.todict(),
        variational = var.todict(),
        minimizer = minimizer.todict(),
        output = output.todict(),
        final = dict(diagnostics=dict(departures='oman')),
        prints = dict(frequency='PT3H')
    )

    yaml.Dumper.ignore_aliases = lambda *args : True
    outfile='test.yml'
    print yaml.dump(data, default_flow_style=False)
    with open('test.yml', 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)

if __name__ == "__main__":
    main()
