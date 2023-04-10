# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
from ewok import JediTask

class forecastMOM6(JediTask):
    def setup(self, config, execs, fix, ic):
        self.command = os.path.join(config['model_path'], "tasks/runForecast.py")