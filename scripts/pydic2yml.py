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
                       Jb={}):
        self.cost_type=cost_type
        self.window_begin=window_begin
        self.window_length=window_length
        self.variables=variables
        self.Jb=Jb

    def todict(self):
        tmpdict = dict(cost_type=self.cost_type, 
                       window_begin=self.window_begin,
                       window_length=self.window_length,
                       variables=self.variables,
                       Jb=self.Jb.todict())
        return tmpdict

class Jb:
    def __init__(self, backgrounds, covariance={}):
        self.backgrounds = backgrounds
        self.covariance = covariance

    def todict(self):
        tmpdict = dict( Background=dict(state=self.backgrounds, 
                        Covariance=self.covariance ) )
        return tmpdict
    
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

#Start of DA window (read from file ...)
yyyy=2015
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

#Date at the middle of the DA window
t_mid=t_start+dt_mid

model_vars=["cicen","hicen","hsnon","tsfcn","qsnon","sicnk","qicnk",
            "socn","tocn","ssh"]
window_begin=t_start.isoformat()+'Z'
window_length='PT'+str(dt.seconds/3600)+'H'


#Setup dict for yml parsing
model=Model()
bkg_mid=BackgroundFile()
bkg_mid=[bkg_mid.todict()]
jb=Jb(backgrounds=bkg_mid)
J=CostFunction(variables=model_vars,Jb=jb)


data = dict(
test_framework_runtime_config = "--log_level=test_suite",
resolution = dict(),
model = model.todict(),
cost_function = J.todict()
)

print yaml.dump(data, default_flow_style=False)
#with open('test.yml', 'w') as outfile:
#    yaml.dump(data, outfile, default_flow_style=False)



