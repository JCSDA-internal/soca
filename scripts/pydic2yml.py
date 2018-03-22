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

class Model:
    def __init__(self, dt="PT1H"):
        self.dt = dt
    
    def todict(self):
        tmpdict = dict( tstep = self.dt )
        return tmpdict



class CostFunction:
    def __init__(self, cost_type='3D-Var',
                       window_begin="2015-01-01T09:00:00Z",
                       window_length="PT6H",
                       variables=["aice","hice"],
                       Jb={},Jo={}):
        self.cost_type=cost_type
        self.window_begin=window_begin
        self.window_length=window_length
        self.variables=variables
        self.Jb=Jb
        self.Jo=Jo

    def todict(self):
        tmpdict = dict(cost_type=self.cost_type, 
                       window_begin=self.window_begin,
                       window_length=self.window_length,
                       variables=self.variables,
                       Jb=self.Jb.todict(),
                       Jo=self.Jo.todict())
        return tmpdict

class Jb:
    def __init__(self, backgrounds, covariance):
        self.backgrounds = backgrounds
        self.covariance = covariance

    def todict(self):
        tmpdict = dict( Background=dict(state=self.backgrounds), 
                        Covariance=self.covariance.todict() )
        return tmpdict

class Jo:
    def __init__(self,obstypes):
        self.obstypes=obstypes

    def todict(self):
        return dict(ObsTypes=obstypes)
    
class ObsType:
    def __init__(self,obstype='SeaIceThickness',
                      obsfileout='Data/test.out',
                      obsfilein='Data/cryosat2-nrt-nh-20171030-v2.0.nc',
                      obsvalue='ObsVal',
                      covariance='diagonal',
                      obserror='ObsErr'):
        self.obstype=obstype
        self.obsfileout=obsfileout
        self.obsfilein=obsfilein
        self.obsvalue=obsvalue
        self.covariance=covariance
        self.obserror=obserror

    def todict(self):
        tmpdict = dict(ObsType=self.obstype,
        ObsData=dict(ObsDataOut=dict(obsfile=self.obsfileout),
                     ObsDataIn=dict(obsfile=self.obsfilein),
                     obsvalue=self.obsvalue),
        Covariance=dict(covariance=self.covariance,
                        obserror=self.obserror))

        return tmpdict        
                 
class Covariance:
    def __init__(self,covariance='SocaError',
                      standard_deviation=0.8,
                      vertical_correlation=0.2,
                      horizontal_length_scale=1e6,
                      date='2015-01-01T12:00:00Z'):
        self.covariance=covariance
        self.standard_deviation=standard_deviation
        self.vertical_correlation=vertical_correlation
        self.horizontal_length_scale=horizontal_length_scale
        self.date=date

    def todict(self):
        return dict(covariance=self.covariance,
                 standard_deviation=self.standard_deviation,
                 vertical_correlation=self.vertical_correlation,
                 horizontal_length_scale=self.horizontal_length_scale,
                 date=self.date)
        
class BackgroundFile:
    def __init__(self,read_from_file=1,
                      basename="./RESTART-01/",
                      ocn_sfc_filename="",
                      ocn_filename="MOM.res.nc",
                      ice_filename="ice_model.res.nc",
                      date="2015-01-01T00:00:00Z"):        
        self.read_from_file=read_from_file
        self.basename=basename
        self.ocn_sfc_filename=ocn_sfc_filename
        self.ocn_filename=ocn_filename
        self.ice_filename=ice_filename
        self.date=date

    def todict(self):
        tmpdict = dict( read_from_file=self.read_from_file,
                      basename=self.basename,
                      ocn_sfc_filename=self.ocn_sfc_filename,
                      ocn_filename=self.ocn_filename,
                      ice_filename=self.ice_filename,
                      date=self.date)
        return tmpdict

class Iteration:
    def __init__(self,resolution,tstep,ninner,tlm='IdTLM'):
        self.resolution=resolution
        self.tstep=tstep
        self.ninner=ninner
        self.tlm=tlm

    def todict(self):
        tmpdict = dict( resolution=self.resolution,
                        linearmodel=dict(version=self.tlm,
                                         tstep=self.tstep),
                        ninner=self.ninner,
                        gradient_norm_reduction=1e-5,
                        test='on' )    
        return tmpdict
    
class Variational:
    def __init__(self,iterations):
        self.iterations=iterations

    def todict(self):
        tmpdict = dict( variational=dict(iteration=self.iterations))
        return tmpdict

class Minimizer:
    def __init__(self,algorithm='DRIPCG'):
        self.algorithm=algorithm

    def todict(self):
        tmpdict = dict(algorithm=self.algorithm)
        return tmpdict

class Output:
    def __init__(self,datadir='Data',exp='3dvar',typeout='an'):
        self.datadir=datadir
        self.exp=exp
        self.typeout=typeout

    def todict(self):
        return dict(datadir=self.datadir,exp=self.exp,type=self.typeout)
    
#Start of DA window (read from file ...)
yyyy=2016
mm=1
dd=1
hh=9
mn=0
ss=0
t_start=datetime.datetime(yyyy,mm,dd,hh,mn,ss)

#DA window length
dt=6 # [hours]
dt_mid=dt/2
dt=datetime.timedelta(hours=dt)
dt_mid=datetime.timedelta(hours=dt_mid)

#Number of outer loop
nouter = 1

#Date at the middle of the DA window
t_mid=t_start+dt_mid
print t_mid

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

print ocnfname, icefname

model_vars=["cicen","hicen","hsnon","tsfcn","qsnon","sicnk","qicnk",
            "socn","tocn","ssh"]
window_begin=t_start.isoformat()+'Z'
window_length='PT'+str(dt.seconds/3600)+'H'

#Setup dict for yml parsing
model=Model()
bkg_error=Covariance()
bkg_mid=BackgroundFile(date=t_mid.isoformat()+'Z',
                       basename='./RESTART/',
                       ocn_filename=ocnfname,
                       ice_filename=icefname)
bkg_mid=[bkg_mid.todict()]
obstypes=ObsType()
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



