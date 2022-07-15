# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import ewok.tasks.getFcInit as generic


class getFcInitMOM6(generic.getFcInit):

    def setup(self, config):
        super().setup(config)

        # Use specific script
        self.command = os.path.join(config['model_path'], "tasks/runGetForecast.py")