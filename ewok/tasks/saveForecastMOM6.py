
# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import ewok.tasks.saveForecast as generic

class saveForecastMOM6(generic.saveForecast):

    def setup(self, config, fc):
        super().setup(config, fc)
        self.command = os.path.join(config['model_path'], "tasks/runSaveForecast.py")