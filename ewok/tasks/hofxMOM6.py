# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
from ewok.tasks.hofx import hofx


class hofxMOM6(hofx):

    def setup(self, config, execs, fix, hcoeffs, bg, obs):
        # overwrite the executable names
        bindir = os.environ.get("JEDI_BUILD")
        execs['hx3dex'] = os.path.join(bindir, 'soca_hofx3d.x')
        execs['hxexec'] = os.path.join(bindir, 'soca_hofx.x')
        super().setup(config, execs, fix, hcoeffs, bg, obs)
