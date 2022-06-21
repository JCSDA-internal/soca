# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
from ewok import Task


class getStaticModelMOM6(Task):

    def setup(self, config, **inputs):
        # localconf = {}
        # localconf['fixdir'] = self.workdir['wdir']
        # localconf['resol'] = config['GEOMETRY']['_resol_name']
        # localconf['nlevs'] = str(config['GEOMETRY']['npz'])

        # ilayout = config['GEOMETRY']['layout']
        # layout = str(ilayout[0]) + "x" + str(ilayout[1])

        # tmplfile = os.path.join(config['model_path'], "templates/fixgeom.yaml")
        # geomtmpl = yamltools.parse_config(tmplfile)
        # geomconf = yamltools.substitute_template_variables(geomtmpl, localconf)
        # self.output['FIXGEOM'] = geomconf

        # tmplfile = os.path.join(config['model_path'], "templates/fixmodel.yaml")
        # modeltmpl = yamltools.parse_config(tmplfile)
        # modelconf = yamltools.substitute_template_variables(modeltmpl, localconf)
        # self.output['FIXMODEL'] = modelconf

        self.RUNTIME_ENV['SOCAREPO'] = os.path.join(os.environ.get("JEDI_SRC"), "soca")
        self.RUNTIME_ENV['SOCA_STATIC_DIR'] = config['soca_static_dir']
        # self.RUNTIME_ENV['WORKDIR'] = self.workdir['wdir']
        # self.RUNTIME_ENV['RESOL'] = config['GEOMETRY']['_resol_name']
        # self.RUNTIME_ENV['NLEVS'] = config['GEOMETRY']['npz']
        # self.RUNTIME_ENV['LAYOUT'] = layout
        self.command = os.path.join(config['model_path'], "tasks/runStaticModelMOM6.sh")