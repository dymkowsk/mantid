# pylint: disable=too-few-public-methods

"""State describing the creation of pixel and wavelength adjustment workspaces for SANS reduction."""

import json
from SANS2.State.SANSStateBase import (SANSStateBase, sans_parameters, StringParameter,
                                       ClassTypeParameter, PositiveFloatParameter, DictParameter)
from SANS2.Common.SANSEnumerations import (RangeStepType, DetectorType, convert_detector_type_to_string)


# ------------------------------------------------
# SANSStateAdjustment
# ------------------------------------------------
class SANSStateWavelengthAndPixelAdjustment(object):
    pass


@sans_parameters
class SANSStateAdjustmentFiles(SANSStateBase):
    pixel_adjustment_file = StringParameter()
    wavelength_adjustment_file = StringParameter()

    def __init__(self):
        super(SANSStateAdjustmentFiles, self).__init__()

    def validate(self):
        is_invalid = {}
        # TODO if a file was specified then make sure that its existence is checked.

        if is_invalid:
            raise ValueError("SANSStateAdjustmentFiles: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))


@sans_parameters
class SANSStateWavelengthAndPixelAdjustmentISIS(SANSStateBase, SANSStateWavelengthAndPixelAdjustment):
    wavelength_low = PositiveFloatParameter()
    wavelength_high = PositiveFloatParameter()
    wavelength_step = PositiveFloatParameter()
    wavelength_step_type = ClassTypeParameter(RangeStepType)

    adjustment_files = DictParameter()

    def __init__(self):
        super(SANSStateWavelengthAndPixelAdjustmentISIS, self).__init__()
        self.adjustment_files = {convert_detector_type_to_string(DetectorType.Lab): SANSStateAdjustmentFiles(),
                                 convert_detector_type_to_string(DetectorType.Hab): SANSStateAdjustmentFiles()}

    def validate(self):
        is_invalid = {}

        if self.wavelength_step_type is None:
            is_invalid.update({"wavelength_step_type": "A wavelength range step type has to be specified."})
        if self.wavelength_step is None:
            is_invalid.update({"wavelength_step": "A wavelength step has to be specified."})
        if self.wavelength_low is None:
            is_invalid.update({"wavelength_low": "A lower wavelength value for rebinning has to be specified."})
        if self.wavelength_high is None:
            is_invalid.update({"wavelength_high": "An high wavelength value for rebinning has to be specified."})
        if self.wavelength_low is not None and self.wavelength_high is not None:
            if self.wavelength_low > self.wavelength_high:
                is_invalid.update({"wavelength_high": "The lower wavelength bound needs to be smaller than the upper "
                                                      "bound. The lower bound is {0} and the upper "
                                                      "is {1}.".format(self.wavelength_low, self.wavelength_high)})

        try:
            self.adjustment_files[convert_detector_type_to_string(DetectorType.Lab)].validate()
            self.adjustment_files[convert_detector_type_to_string(DetectorType.Hab)].validate()
        except ValueError as e:
            is_invalid.update({"adjustment_files": str(e)})

        if is_invalid:
            raise ValueError("SANSStateWavelengthAndPixelAdjustmentISIS: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))


# -----------------------------------------------
# SANSStateNormalizeMonitor setup for other facilities/techniques/scenarios.
# Needs to derive from SANSStateNormalizeMonitor and SANSStateBase and fulfill its contract.
# -----------------------------------------------