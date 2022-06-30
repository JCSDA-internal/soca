# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
from ewok import Task


class getStaticModelMOM6(Task):

    def setup(self, config, **inputs):

        self.RUNTIME_ENV['SOCAREPO'] = os.path.join(os.environ.get("JEDI_SRC"), "soca")
        self.RUNTIME_ENV['SOCA_STATIC_DIR'] = config['soca_static_dir']
        self.RUNTIME_ENV['RESOLUTION'] = config['GEOMETRY']['_resol_name']

        self.command = os.path.join(config['model_path'], "tasks/runStaticModelMOM6.sh")