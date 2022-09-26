
# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import ewok.tasks.saveAnalysis as generic

class saveAnalysisMOM6(generic.saveAnalysis):

    def setup(self, config, an):
        super().setup(config, an)
        self.command = os.path.join(config['model_path'], "tasks/runSaveAnalysis.py")