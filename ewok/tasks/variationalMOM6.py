# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
from ewok.tasks.variational import variational

class variationalMOM6(variational):
    def setup(self, config, execs, fix, hcoeffs, statb, bg, bgens, obs):
        # overwrite the executable names
        build = os.environ.get("JEDI_BUILD")
        bindir = os.path.join(build, 'bin')
        execs['anexec'] = os.path.join(bindir, 'soca_var.x')
        super().setup(config, execs, fix, hcoeffs, statb, bg, bgens, obs)
