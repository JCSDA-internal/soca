
# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import ewok
from ewok.tasks import GenericModel

from forecastMOM6 import forecastMOM6
from getBackgroundMOM6 import getBackgroundMOM6
from getFCInitMOM6 import getFcInitMOM6
from getInitialConditionsMOM6 import getInitialConditionsMOM6
from getStaticModelMOM6 import getStaticModelMOM6
from hofxMOM6 import hofxMOM6
from saveAnalysisMOM6 import saveAnalysisMOM6
from saveForecastMOM6 import saveForecastMOM6
from variationalMOM6 import variationalMOM6

class ModelTasks(GenericModel.ModelTasks):
    def __init__(self):
        super().__init__()

        # self.getFcInit = getFcInitMOM6
        self.createPlots = ewok.createPlots
        self.forecast = forecastMOM6
        self.getBackground = getBackgroundMOM6
        self.getInitialConditions = getInitialConditionsMOM6
        self.getStaticModel = getStaticModelMOM6
        self.hofx = hofxMOM6
        self.saveAnalysis = saveAnalysisMOM6
        self.saveForecast = saveForecastMOM6
        self.savePlots = ewok.savePlots
        self.variational = variationalMOM6