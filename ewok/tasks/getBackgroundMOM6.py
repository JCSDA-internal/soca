
# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.getBackground as generic


class getBackgroundMOM6(generic.getBackground):

    def setup(self, config, fc):
        super().setup(config, fc)

        # Use MOM6 specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetForecast.py")