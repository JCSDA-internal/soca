
# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

from ewok.tasks import GenericModel

from hofxMOM6 import hofxMOM6
from getBackgroundMOM6 import getBackgroundMOM6
from getFCInitMOM6 import getFcInitMOM6
from getInitialConditionsMOM6 import getInitialConditionsMOM6
from getStaticModelMOM6 import getStaticModelMOM6

class ModelTasks(GenericModel.ModelTasks):
    def __init__(self):
        super().__init__()

        # self.getBackground = getBackgroundMOM6
        # self.getFcInit = getFcInitMOM6
        self.getInitialConditions = getInitialConditionsMOM6
        self.hofx = hofxMOM6
        self.getStaticModel = getStaticModelMOM6