# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import ewok.tasks.getInitialConditions as generic


class getInitialConditionsMOM6(generic.getInitialConditions):

    def setup(self, config):
        super().setup(config)

        # self.walltime = '00:05:00'

        # Use specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetForecast.py")