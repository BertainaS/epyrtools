# fair_converter.py
import sys
import os
import json
import csv
from pathlib import Path
import numpy as np
import warnings
import h5py  # Required for HDF5 output: pip install h5py
from typing import Dict, Any, Tuple, Optional, Union, List
from datetime import datetime # For timestamp

# Assuming eprload is in the same directory or installable
try:
    from .eprload import eprload
except ImportError:
    print("Error: Could not import 'eprload'. Make sure eprload.py is in the same directory or accessible in the Python path.", file=sys.stderr)
    sys.exit(1)

BRUKER_PARAM_MAP = {
    # ----- DSC Section (Data Set Codes) -----
    "DSRC": {"fair_name": "data_source_type", "unit": None, "description": "Type of dataset (e.g., EXP for experimental data)."},
    "BSEQ": {"fair_name": "byte_sequence_order", "unit": None, "description": "Byte sequence order (e.g., BIG for big-endian)."},
    "IKKF": {"fair_name": "imaginary_data_kind_format", "unit": None, "description": "Format of imaginary data component (e.g., REAL, CMPLX, NODATA)."},
    "XTYP": {"fair_name": "x_axis_type", "unit": None, "description": "Type of X-axis data (e.g., IDX for indexed, MONO for monotonic)."},
    "YTYP": {"fair_name": "y_axis_type", "unit": None, "description": "Type of Y-axis data (e.g., IDX for indexed, NODATA for 1D spectra)."},
    "ZTYP": {"fair_name": "z_axis_type", "unit": None, "description": "Type of Z-axis data (e.g., NODATA for 1D/2D spectra)."},
    "IRFMT": {"fair_name": "intensity_data_format", "unit": None, "description": "Format of real part of intensity data (e.g., D for double precision float)."},
    "IIFMT": {"fair_name": "imaginary_intensity_data_format", "unit": None, "description": "Format of imaginary part of intensity data (e.g., D for double precision float)."}, # New
    "XPTS": {"fair_name": "number_of_points_x_axis", "unit": "points", "description": "Number of data points along the X-axis."},
    "XMIN": {"fair_name": "x_axis_minimum", "unit": "refer to XUNI", "description": "Minimum value of the X-axis. Unit specified by XUNI."},
    "XWID": {"fair_name": "x_axis_width", "unit": "refer to XUNI", "description": "Width (range) of the X-axis. Unit specified by XUNI."},
    "YPTS": {"fair_name": "number_of_points_y_axis", "unit": "points", "description": "Number of data points along the Y-axis (for 2D experiments)."}, # New
    "YMIN": {"fair_name": "y_axis_minimum", "unit": "refer to YUNI", "description": "Minimum value of the Y-axis (for 2D experiments). Unit specified by YUNI."}, # New
    "YWID": {"fair_name": "y_axis_width", "unit": "refer to YUNI", "description": "Width (range) of the Y-axis (for 2D experiments). Unit specified by YUNI."}, # New
    "TITL": {"fair_name": "title", "unit": None, "description": "Title of the experiment or dataset."},
    "IRNAM": {"fair_name": "intensity_axis_name_real", "unit": None, "description": "Name of the real part of the intensity (ordinate) axis."},
    "IINAM": {"fair_name": "intensity_axis_name_imaginary", "unit": None, "description": "Name of the imaginary part of the intensity (ordinate) axis."}, # New
    "XNAM": {"fair_name": "x_axis_name", "unit": None, "description": "Name of the X-axis (abscissa)."},
    "YNAM": {"fair_name": "y_axis_name", "unit": None, "description": "Name of the Y-axis (second dimension)."}, # New
    "IRUNI": {"fair_name": "intensity_axis_unit_real", "unit": None, "description": "Unit of the real part of the intensity (ordinate) axis."},
    "IIUNI": {"fair_name": "intensity_axis_unit_imaginary", "unit": None, "description": "Unit of the imaginary part of the intensity (ordinate) axis."}, # New
    "XUNI": {"fair_name": "x_axis_unit", "unit": None, "description": "Unit of the X-axis (abscissa)."},
    "YUNI": {"fair_name": "y_axis_unit", "unit": None, "description": "Unit of the Y-axis (second dimension)."}, # New

    # SPL Section (Standard Parameter Layer)
    "OPER": {"fair_name": "operator_name", "unit": None, "description": "Name or identifier of the operator."},
    "DATE": {"fair_name": "acquisition_date", "unit": None, "description": "Date of the experiment."},
    "TIME": {"fair_name": "acquisition_time", "unit": None, "description": "Time of the experiment."},
    "CMNT": {"fair_name": "comment", "unit": None, "description": "User-defined comment for the experiment."},
    "SAMP": {"fair_name": "sample_identifier", "unit": None, "description": "Sample name or identifier."},
    "SFOR": {"fair_name": "sample_form", "unit": None, "description": "Description of the sample's form or preparation."},
    "STAG": {"fair_name": "sample_stage_temperature", "unit": "K", "description": "Sample stage or temperature controller setting. Assumed Kelvin, verify units."}, # Often temperature
    "EXPT": {"fair_name": "experiment_type", "unit": None, "description": "Type of experiment (e.g., CW for Continuous Wave, PLS for Pulsed)."},
    "OXS1": {"fair_name": "ordinate_axis_1_source", "unit": None, "description": "Source/type for the first ordinate axis (e.g., IADC - Imaginary ADC signal, TADC - Transient ADC)."},
    "AXS1": {"fair_name": "abscissa_axis_1_type", "unit": None, "description": "Type/parameter for the first abscissa axis (e.g., B0VL - Magnetic Field)."},
    "AXS2": {"fair_name": "abscissa_axis_2_type", "unit": None, "description": "Type/parameter for the second abscissa axis (if used)."},
    "AXS3": {"fair_name": "abscissa_axis_3_type", "unit": None, "description": "Type/parameter for the third abscissa axis (if used)."},
    "A1CT": {"fair_name": "abscissa_1_center", "unit": None, "description": "Center value for the first abscissa axis. Units depend on AXS1 (e.g., Tesla if B0VL)."},
    "A1SW": {"fair_name": "abscissa_1_sweep_width", "unit": None, "description": "Sweep width for the first abscissa axis. Units depend on AXS1 (e.g., Tesla if B0VL)."},
    "MWFQ": {"fair_name": "microwave_frequency", "unit": "Hz", "description": "Microwave frequency."},
    "MWPW": {"fair_name": "microwave_power", "unit": "W", "description": "Microwave power."}, # Changed from mW to W for SPL; DSL might use mW
    "AVGS": {"fair_name": "number_of_averages_set", "unit": "scans", "description": "Number of averages/scans set for accumulation."},
    "RESO": {"fair_name": "resonator_identifier", "unit": None, "description": "Identifier or type of the resonator used."},
    "SPTP": {"fair_name": "sweep_time_per_point", "unit": "s", "description": "Sweep time per data point (dwell time)."},
    "RCAG": {"fair_name": "receiver_gain_spl", "unit": "dB", "description": "Receiver gain value from Standard Parameter Layer."},
    "RCHM": {"fair_name": "receiver_harmonic_spl", "unit": None, "description": "Receiver harmonic setting from Standard Parameter Layer."},
    "B0MA": {"fair_name": "modulation_amplitude_tesla", "unit": "T", "description": "Magnetic field modulation amplitude in Tesla."},
    "B0MF": {"fair_name": "modulation_frequency_spl", "unit": "Hz", "description": "Magnetic field modulation frequency in Hz from Standard Parameter Layer."},
    "RCPH": {"fair_name": "receiver_phase", "unit": "deg", "description": "Receiver phase."},
    "RCOF": {"fair_name": "receiver_offset_spl", "unit": "%", "description": "Receiver offset as a percentage from Standard Parameter Layer."},
    "A1RS": {"fair_name": "abscissa_1_resolution_points", "unit": "points", "description": "Resolution (number of points) for the first abscissa axis."},
    "RCTC": {"fair_name": "receiver_time_constant_spl", "unit": "s", "description": "Receiver time constant in seconds from Standard Parameter Layer."},
    "B0VL": {"fair_name": "magnetic_field_value", "unit": "T", "description": "Magnetic field value."},  # New parameter from SPL

    # DSL Section - Parameters are prefixed with their DVC block for clarity
    # .DVC acqStart
    "acqStart.comment": {"fair_name": "acquisition_start_comment", "unit": None, "description": "Comment related to the start of acquisition."},

    # .DVC cwBridge
    "AcqFineTuning": {"fair_name": "cw_bridge_acquisition_fine_tuning", "unit": None, "description": "Continuous Wave bridge fine tuning mode during acquisition."},
    "AcqScanFTuning": {"fair_name": "cw_bridge_acquisition_scan_fine_tuning", "unit": None, "description": "Continuous Wave bridge scan fine tuning status (On/Off)."},
    "AcqSliceFTuning": {"fair_name": "cw_bridge_acquisition_slice_fine_tuning", "unit": None, "description": "Continuous Wave bridge slice fine tuning status (On/Off)."},
    "BridgeCalib": {"fair_name": "cw_bridge_calibration_value", "unit": None, "description": "Bridge calibration value."}, # Unit might be dB or similar, needs checking
    "Power": {"fair_name": "cw_bridge_power_output", "unit": "mW", "description": "Microwave power setting at the CW bridge output."},
    "PowerAtten": {"fair_name": "cw_bridge_power_attenuation", "unit": "dB", "description": "Microwave power attenuation setting at the CW bridge."},

    # .DVC endor (Electron Nuclear Double Resonance)
    "EIEENDORFreq": {"fair_name": "endor_eie_frequency", "unit": "MHz/kG", "description": "ENDOR EIE (Electron Irradiation Effect) frequency per kG."},
    "EIEIsotope": {"fair_name": "endor_eie_isotope", "unit": None, "description": "Isotope used for ENDOR EIE measurement."},
    "RFSweepDir": {"fair_name": "endor_rf_sweep_direction", "unit": None, "description": "Direction of RF sweep for ENDOR (e.g., Same, Opposite)."},
    "EIEStaticField": {"fair_name": "endor_eie_static_field", "unit": "G", "description": "Static magnetic field for ENDOR EIE."},
    "EIEStaticRF": {"fair_name": "endor_eie_static_rf_frequency", "unit": "MHz", "description": "Static RF frequency for ENDOR EIE."},
    "ENDORType": {"fair_name": "endor_type", "unit": None, "description": "Type of ENDOR experiment (e.g., EIF)."},
    "RF1Atten": {"fair_name": "endor_rf1_attenuation", "unit": "dB", "description": "Attenuation for RF channel 1."},
    "RF1FreqPos": {"fair_name": "endor_rf1_frequency_position", "unit": "MHz", "description": "Frequency position for RF channel 1."},
    "RF1StartFreq": {"fair_name": "endor_rf1_start_frequency", "unit": "MHz", "description": "Start frequency for RF channel 1 sweep."},
    "RF1SweepWidth": {"fair_name": "endor_rf1_sweep_width", "unit": "MHz", "description": "Sweep width for RF channel 1."},
    "RF2Atten": {"fair_name": "endor_rf2_attenuation", "unit": "dB", "description": "Attenuation for RF channel 2."},
    "RF2FreqPos": {"fair_name": "endor_rf2_frequency_position", "unit": "MHz", "description": "Frequency position for RF channel 2."},
    "RF2StartFreq": {"fair_name": "endor_rf2_start_frequency", "unit": "MHz", "description": "Start frequency for RF channel 2 sweep."},
    "RF2SweepWidth": {"fair_name": "endor_rf2_sweep_width", "unit": "MHz", "description": "Sweep width for RF channel 2."},
    "RFSrcMixing": {"fair_name": "endor_rf_source_mixing", "unit": None, "description": "RF source mixing mode (e.g., Add)."},
    "SumAtten": {"fair_name": "endor_sum_attenuation", "unit": "dB", "description": "Summation attenuator setting for ENDOR."},
    "SumAttenStart": {"fair_name": "endor_sum_attenuation_start", "unit": "dB", "description": "Start value for sum attenuation sweep."},
    "SumAttenWidth": {"fair_name": "endor_sum_attenuation_width", "unit": "dB", "description": "Width for sum attenuation sweep."},

    # .DVC fieldCtrl (assuming some overlap with previous example, adding new ones)
    "AllegroMode": {"fair_name": "field_controller_allegro_mode", "unit": None, "description": "Allegro mode status for field controller."},
    "CenterField": {"fair_name": "field_controller_center_field", "unit": "G", "description": "Center magnetic field setting for the field controller."},
    "Delay": {"fair_name": "field_controller_delay", "unit": "s", "description": "Delay after setting field before measurement by field controller."},
    "FieldFlyback": {"fair_name": "field_controller_flyback_enabled", "unit": None, "description": "Field flyback status (On/Off) for field controller."},
    "FieldPosition": {"fair_name": "field_controller_current_field_position", "unit": "G", "description": "Current magnetic field position reported by the field controller."}, # New
    "FieldWait": {"fair_name": "field_controller_wait_condition", "unit": None, "description": "Condition for waiting after field set by field controller (e.g., 'Wait LED off')."},
    "GFactor": {"fair_name": "field_controller_g_factor_reference", "unit": None, "description": "g-factor setting used for field-frequency calculations by field controller."},
    "MeasuringHall": {"fair_name": "field_controller_measuring_hall_sensor", "unit": None, "description": "Indicates if the Hall sensor is actively measuring."},
    "SetToSampleG": {"fair_name": "field_controller_set_to_sample_g_factor", "unit": None, "description": "Indicates if field controller uses a sample-specific g-factor."},
    "StaticFieldMon": {"fair_name": "field_controller_static_field_monitor", "unit": "G", "description": "Monitored static magnetic field value."},
    "SweepDirection": {"fair_name": "field_controller_sweep_direction", "unit": None, "description": "Direction of magnetic field sweep by field controller (e.g., Up, Down)."},
    "SweepWidth": {"fair_name": "field_controller_sweep_width", "unit": "G", "description": "Magnetic field sweep width setting for the field controller."},
    "WidthTM": {"fair_name": "field_controller_width_tm", "unit": "G", "description": "Teslameter field range or similar, context needed."}, # Unclear, assuming related to field measurement range

    # .DVC freqCounter
    "FrequencyMon": {"fair_name": "frequency_counter_monitored_frequency", "unit": "GHz", "description": "Microwave frequency as monitored by the frequency counter."},
    "QMonitBridge": {"fair_name": "frequency_counter_q_monitor_bridge_enabled", "unit": None, "description": "Q-factor monitoring bridge status via frequency counter."},

    # .DVC ftBridge (Assuming Fourier Transform Bridge specific parameters)
    "Attenuation": {"fair_name": "ft_bridge_attenuation", "unit": "dB", "description": "Attenuation setting for the FT bridge."},
    "ELDORAtt": {"fair_name": "ft_bridge_eldor_attenuation", "unit": "dB", "description": "ELDOR attenuation setting within the FT bridge."},
    "FrequencyA": {"fair_name": "ft_bridge_frequency_a", "unit": "GHz", "description": "Frequency setting for channel A of the FT bridge."},
    "VideoBW": {"fair_name": "ft_bridge_video_bandwidth", "unit": "MHz", "description": "Video bandwidth setting for the FT bridge."},
    "VideoGain": {"fair_name": "ft_bridge_video_gain", "unit": "dB", "description": "Video gain setting for the FT bridge."},

    # .DVC ftEpr (Fourier Transform EPR specific parameters)
    "AWGPhaseShift": {"fair_name": "ftepr_awg_phase_shift", "unit": "deg", "description": "Phase shift applied by the Arbitrary Waveform Generator."},
    "AWGPrg": {"fair_name": "ftepr_awg_program", "unit": None, "description": "Arbitrary Waveform Generator program string."},
    "AutoTimeOut": {"fair_name": "ftepr_auto_timeout_enabled", "unit": None, "description": "Automatic timeout status for FT-EPR."},
    "AveragesPerScan": {"fair_name": "ftepr_averages_per_scan", "unit": None, "description": "Number of averages per scan in FT-EPR."},
    "ELDORFreqStart": {"fair_name": "ftepr_eldor_frequency_start", "unit": "GHz", "description": "Start frequency for ELDOR in FT-EPR."},
    "ELDORFreqWidth": {"fair_name": "ftepr_eldor_frequency_width", "unit": "GHz", "description": "Frequency width for ELDOR in FT-EPR."},
    "FTEzAWGELDORa": {"fair_name": "ftepr_awg_eldor_amplitude", "unit": "%", "description": "Amplitude of the ELDOR pulse from AWG."},
    "FTEzAWGELDORf": {"fair_name": "ftepr_awg_eldor_frequency_offset", "unit": "MHz", "description": "Frequency offset for the ELDOR pulse from AWG."},
    "FTEzAWGELDORw": {"fair_name": "ftepr_awg_eldor_pulse_width", "unit": "MHz", "description": "Pulse width for the ELDOR pulse from AWG."}, # MHz might be a typo for ns or us; needs verification
    "FieldIsStatic": {"fair_name": "ftepr_field_is_static", "unit": None, "description": "Indicates if the magnetic field is static during the experiment."},
    "GradIntPulse": {"fair_name": "ftepr_gradient_integration_pulse", "unit": None, "description": "Status of gradient integration pulse."},
    "GrdEnable": {"fair_name": "ftepr_gradient_enable", "unit": None, "description": "Status of gradient enable."},
    "LastXAxis": {"fair_name": "ftepr_last_x_axis_type", "unit": None, "description": "Type of the last X-axis used."},
    "LastYAxis": {"fair_name": "ftepr_last_y_axis_type", "unit": None, "description": "Type of the last Y-axis used."},
    "MDigPrg": {"fair_name": "ftepr_mdig_program", "unit": None, "description": "Magnet Digital Controller program string."},
    "MMWaveLOFreq": {"fair_name": "ftepr_mmwave_local_oscillator_frequency", "unit": "GHz", "description": "Millimeter-wave local oscillator frequency."},
    "MicroImgPrg": {"fair_name": "ftepr_micro_imaging_program", "unit": None, "description": "Micro-imaging program string."},
    "OnlyAWGChans": {"fair_name": "ftepr_only_awg_channels", "unit": None, "description": "Indicates if only AWG channels are used."},
    "PCycleAllowed": {"fair_name": "ftepr_phase_cycling_allowed", "unit": None, "description": "Indicates if phase cycling is allowed."},
    "PCycleOn": {"fair_name": "ftepr_phase_cycling_enabled", "unit": None, "description": "Indicates if phase cycling is enabled."},
    "PPExtTrg": {"fair_name": "ftepr_pulse_program_external_trigger", "unit": None, "description": "Status of external trigger for pulse program."},
    "PPExtTrgSlope": {"fair_name": "ftepr_pulse_program_external_trigger_slope", "unit": None, "description": "Slope of the external trigger for pulse program."},
    "PlsSPELEXPSlct": {"fair_name": "ftepr_pulse_program_experiment_select", "unit": None, "description": "Selected pulse program/experiment."},
    "PlsSPELGlbTxt": {"fair_name": "ftepr_pulse_program_global_text", "unit": None, "description": "Global text or comments for the pulse program."},
    "PlsSPELLISTSlct": {"fair_name": "ftepr_pulse_program_list_select", "unit": None, "description": "Selected list for the pulse program."},
    "PlsSPELPhPrgEx": {"fair_name": "ftepr_pulse_program_phase_program_extended", "unit": None, "description": "Extended phase program for the pulse sequence."},
    "PlsSPELPrg": {"fair_name": "ftepr_pulse_program_name", "unit": None, "description": "Name of the pulse program file."},
    "PlsSPELShpTxt": {"fair_name": "ftepr_pulse_shape_text", "unit": None, "description": "Text defining pulse shapes."},
    "Psd1": {"fair_name": "ftepr_pulse_sequence_definition_1", "unit": None, "description": "Pulse sequence definition part 1."}, # These Psd parameters likely define pulse sequence elements.
    "Psd10": {"fair_name": "ftepr_pulse_sequence_definition_10", "unit": None, "description": "Pulse sequence definition part 10."},
    "Psd11": {"fair_name": "ftepr_pulse_sequence_definition_11", "unit": None, "description": "Pulse sequence definition part 11."},
    "PlsSPELPrgTxt": {"fair_name": "ftepr_pulse_program_text", "unit": None, "description": "Full text of the pulse program."},
    "TD_ADC_IntegEnd": {"fair_name": "ftepr_transient_adc_integration_end", "unit": "ns", "description": "End time for ADC integration in transient mode."},
    "TD_ADC_IntegStart": {"fair_name": "ftepr_transient_adc_integration_start", "unit": "ns", "description": "Start time for ADC integration in transient mode."},
    "TD_ArrayValues1": {"fair_name": "ftepr_transient_array_values_1", "unit": None, "description": "Array values for transient experiment, dimension 1."},
    "TD_Background": {"fair_name": "ftepr_transient_background_subtraction", "unit": None, "description": "Background subtraction setting for transient data."},
    "TD_CaptureLen": {"fair_name": "ftepr_transient_capture_length", "unit": "samples", "description": "Length of the captured transient data."},
    "TD_DelayToSampling": {"fair_name": "ftepr_transient_delay_to_sampling", "unit": "ns", "description": "Delay before starting data sampling in transient mode."},
    "TD_DirectAxis": {"fair_name": "ftepr_transient_direct_axis", "unit": None, "description": "Specifies the direct acquisition dimension in transient experiments."},
    "TD_EnableFIR": {"fair_name": "ftepr_transient_fir_filter_enabled", "unit": None, "description": "Status of FIR (Finite Impulse Response) filter."},
    "TD_EnableScaling": {"fair_name": "ftepr_transient_scaling_enabled", "unit": None, "description": "Status of data scaling for transient data."},
    "TD_NumPoints": {"fair_name": "ftepr_transient_number_of_points", "unit": "points", "description": "Number of points acquired in the transient."},
    "TD_PostTrigTime": {"fair_name": "ftepr_transient_post_trigger_time", "unit": "ns", "description": "Post-trigger time for transient acquisition."},
    "TD_ScaleFactor": {"fair_name": "ftepr_transient_scaling_factor", "unit": None, "description": "Scaling factor applied to transient data."},
    "TD_TimeOffset": {"fair_name": "ftepr_transient_time_offset", "unit": "ns", "description": "Time offset for transient acquisition."},
    "TD_TriggerPosition": {"fair_name": "ftepr_transient_trigger_position", "unit": "%", "description": "Trigger position as a percentage of the acquisition window."},

    # .DVC recorder (some might overlap with previous example)
    "BaselineCorr": {"fair_name": "recorder_baseline_correction", "unit": None, "description": "Recorder baseline correction setting."},
    "NbScansAcc": {"fair_name": "recorder_accumulated_scans", "unit": None, "description": "Number of scans accumulated by the recorder."},
    "NbScansDone": {"fair_name": "recorder_scans_completed", "unit": None, "description": "Number of scans completed by the recorder."},
    "NbScansToDo": {"fair_name": "recorder_scans_to_do", "unit": None, "description": "Number of scans remaining to be performed."},
    "ReplaceMode": {"fair_name": "recorder_replace_mode", "unit": None, "description": "Recorder data replace mode (On/Off)."},
    "SmoothMode": {"fair_name": "recorder_smooth_mode", "unit": None, "description": "Smoothing mode applied by the recorder (e.g., Manual)."},
    "SmoothPoints": {"fair_name": "recorder_smooth_points", "unit": None, "description": "Number of points used for smoothing."},

    # .DVC signalChannel (some might overlap with previous example)
    "AFCTrap": {"fair_name": "signal_channel_afc_trap", "unit": None, "description": "Automatic Frequency Control (AFC) trap status."},
    "Calibrated": {"fair_name": "signal_channel_calibrated", "unit": None, "description": "Indicates if the signal channel is calibrated."},
    "ConvTime": {"fair_name": "signal_channel_conversion_time", "unit": "ms", "description": "Analog-to-digital conversion time."},
    "DModAFCTrap": {"fair_name": "signal_channel_double_modulation_afc_trap", "unit": None, "description": "Double modulation AFC trap status."},
    "DModAmp": {"fair_name": "signal_channel_double_modulation_amplitude", "unit": "G", "description": "Double modulation amplitude."},
    "DModCalibrated": {"fair_name": "signal_channel_double_modulation_calibrated", "unit": None, "description": "Calibration status for double modulation."},
    "DModDetectSCT": {"fair_name": "signal_channel_double_modulation_detection_sct", "unit": None, "description": "Double modulation detection SCT (Signal Channel Transducer) setting."},
    "DModEliDelay": {"fair_name": "signal_channel_double_modulation_elimination_delay", "unit": "us", "description": "Elimination delay for double modulation."},
    "DModExtLockIn": {"fair_name": "signal_channel_double_modulation_external_lockin", "unit": None, "description": "External lock-in status for double modulation."},
    "DModExtTrigger": {"fair_name": "signal_channel_double_modulation_external_trigger", "unit": None, "description": "External trigger status for double modulation."},
    "DModFieldMod": {"fair_name": "signal_channel_double_modulation_field_modulation", "unit": None, "description": "Field modulation type for double modulation."},
    "DModGain": {"fair_name": "signal_channel_double_modulation_gain", "unit": "dB", "description": "Gain for double modulation."},
    "DModHighPass": {"fair_name": "signal_channel_double_modulation_high_pass_filter", "unit": None, "description": "High-pass filter status for double modulation."},
    "DModIntegrator": {"fair_name": "signal_channel_double_modulation_integrator", "unit": None, "description": "Integrator status for double modulation."},
    "DModModOutput": {"fair_name": "signal_channel_double_modulation_output", "unit": None, "description": "Modulation output source for double modulation (e.g., Internal)."},
    "DModSignalInput": {"fair_name": "signal_channel_double_modulation_signal_input", "unit": None, "description": "Signal input source for double modulation (e.g., Internal)."},
    "DModTimeConst": {"fair_name": "signal_channel_double_modulation_time_constant", "unit": "ms", "description": "Time constant for double modulation."},
    "DoubleMode": {"fair_name": "signal_channel_double_modulation_mode_enabled", "unit": None, "description": "Double modulation mode status."},
    "EliDelay": {"fair_name": "signal_channel_elimination_delay", "unit": "us", "description": "Elimination delay in the signal channel."},
    "ExtLockIn": {"fair_name": "signal_channel_external_lockin_enabled", "unit": None, "description": "External lock-in status."},
    "ExtTrigger": {"fair_name": "signal_channel_external_trigger_enabled", "unit": None, "description": "External trigger status."},
    "Gain": {"fair_name": "signal_channel_gain", "unit": "dB", "description": "Gain of the signal channel."},
    "Harmonic": {"fair_name": "signal_channel_harmonic", "unit": None, "description": "Harmonic number for detection."},
    "HighPass": {"fair_name": "signal_channel_high_pass_filter_enabled", "unit": None, "description": "High-pass filter status."},
    "Integrator": {"fair_name": "signal_channel_integrator_enabled", "unit": None, "description": "Integrator status."},
    "ModAmp": {"fair_name": "signal_channel_modulation_amplitude", "unit": "G", "description": "Modulation amplitude."},
    "ModFreq": {"fair_name": "signal_channel_modulation_frequency", "unit": "kHz", "description": "Modulation frequency."},
    "ModInput": {"fair_name": "signal_channel_modulation_input_source", "unit": None, "description": "Modulation input source (e.g., Internal)."},
    "ModOutput": {"fair_name": "signal_channel_modulation_output_source", "unit": None, "description": "Modulation output source (e.g., Internal)."},
    "ModPhase": {"fair_name": "signal_channel_modulation_phase", "unit": "deg", "description": "Modulation phase."},
    "Offset": {"fair_name": "signal_channel_offset", "unit": "%", "description": "Signal offset percentage."},
    "QuadMode": {"fair_name": "signal_channel_quadrature_mode_enabled", "unit": None, "description": "Quadrature detection mode status."},
    "QuadPhase": {"fair_name": "signal_channel_quadrature_phase", "unit": "deg", "description": "Quadrature phase."},
    "Resolution": {"fair_name": "signal_channel_resolution", "unit": "points", "description": "Signal acquisition resolution (number of points)."}, # Note: A1RS is often the primary source for points in X
    "Resonator": {"fair_name": "signal_channel_resonator_selection", "unit": None, "description": "Selected resonator ID or type."},
    "SctNorm": {"fair_name": "signal_channel_sct_normalization", "unit": None, "description": "Signal Channel Transducer normalization status."},
    "SctRevision": {"fair_name": "signal_channel_sct_revision", "unit": None, "description": "Signal Channel Transducer revision."},
    "SignalInput": {"fair_name": "signal_channel_input_source", "unit": None, "description": "Signal input source (e.g., Internal)."},
    "SweepTime": {"fair_name": "signal_channel_sweep_time_actual", "unit": "s", "description": "Actual sweep time for the signal channel."},
    "TimeConst": {"fair_name": "signal_channel_time_constant", "unit": "ms", "description": "Time constant for the signal channel."},
    "TuneCaps": {"fair_name": "signal_channel_tuning_capacitors_value", "unit": None, "description": "Tuning capacitors value or setting."},

    # .DVC transRec (Transient Recorder)
    "transRec.Delay": {"fair_name": "transient_recorder_delay", "unit": "ns", "description": "Delay before transient recording."}, # Assuming ns, needs verification
    "transRec.Length": {"fair_name": "transient_recorder_length", "unit": "ns", "description": "Length of the transient recording."}, # Assuming ns, needs verification
    "transRec.Position": {"fair_name": "transient_recorder_position", "unit": "%", "description": "Position of the transient recording window relative to the trigger (percentage)."},

    # DSC Parameters with Units (Examples - needs verification of actual parameter names and context)
    "XMIN": {"fair_name": "x_axis_minimum", "unit": "G", "description": "Minimum value of the X-axis (unit defined by XUNI)."},
    "XWID": {"fair_name": "x_axis_width", "unit": "G", "description": "Width (range) of the X-axis (unit defined by XUNI)."},
    "YMIN": {"fair_name": "y_axis_minimum", "unit": "MHz", "description": "Minimum value of the Y-axis (unit defined by YUNI)."}, # Assuming Y is frequency here
    "YWID": {"fair_name": "y_axis_width", "unit": "MHz", "description": "Width (range) of the Y-axis (unit defined by YUNI)."}, # Assuming Y is frequency here

    # SPL Parameters with Units (Examples - needs verification of actual parameter names and context)
    "MWFQ": {"fair_name": "microwave_frequency", "unit": "GHz", "description": "Microwave frequency."}, # Common notation is GHz
    "MWPW": {"fair_name": "microwave_power", "unit": "mW", "description": "Microwave power."},
    "B0MA": {"fair_name": "modulation_amplitude", "unit": "G", "description": "Magnetic field modulation amplitude in Gauss."},
    "RCTC": {"fair_name": "receiver_time_constant", "unit": "ms", "description": "Receiver time constant in milliseconds."}, # Often in ms

    # ... (continue adding more based on typical Bruker DSC/DSC files)
  # --- ESP Parameters ---
  # --- General & File ---
  "DOS  Format": { # Note: Two spaces after DOS if that's in your file
    "fair_name": "file_format_version",
    "unit": None,
    "description": "Specifies the DOS text file format version for the parameter file."
  },
  "JON": {
    "fair_name": "instrument_manufacturer_jcamp",
    "unit": None,
    "description": "JCAMP-DX originating system or instrument manufacturer."
  },
  "JDA": {
    "fair_name": "acquisition_date_jcamp",
    "unit": None, # String representation of date
    "description": "JCAMP-DX date of data acquisition (DD/Mon/YYYY)."
  },
  "JTM": {
    "fair_name": "acquisition_time_jcamp",
    "unit": None, # String representation of time
    "description": "JCAMP-DX time of data acquisition (HH:MM)."
  },

  # --- Data Dimensions & Points ---
  "ANZ": {
    "fair_name": "total_number_of_points", # For 2D, this is typically points_dim1 * points_dim2
    "unit": None, # Typically an integer count
    "description": "Total number of data points in the dataset (e.g., X_points * Y_points for 2D)."
  },
  "MIN": {
    "fair_name": "data_axis_minimum_value_defined", # For 2D, this applies to the linearized data
    "unit": "G", # Or could be unitless depending on context, example showed G
    "description": "Defined minimum value of the data axis. For 2D, this may relate to the full range of the ADC for the linearized 2D data block or a display default. Actual swept range per dimension is derived from other parameters."
  },
  "MAX": {
    "fair_name": "data_axis_maximum_value_defined", # For 2D, this applies to the linearized data
    "unit": "G", # Or could be unitless, example showed G
    "description": "Defined maximum value of the data axis. For 2D, this may relate to the full range of the ADC for the linearized 2D data block or a display default. Actual swept range per dimension is derived from other parameters."
  },
  "JSS": { # Interpretation might vary for 2D
    "fair_name": "jcamp_spectral_data_structure_parameter",
    "unit": None,
    "description": "JCAMP-DX parameter related to spectral data structure. For 1D, spectrum status. For 2D, might be points in first dimension (e.g. F2 size) or data block size."
  },
  "SSX": {
    "fair_name": "number_of_points_x_dimension",
    "unit": None,
    "description": "Number of data points in the X-dimension (often the directly acquired dimension or F2)."
  },
  "SSY": {
    "fair_name": "number_of_points_y_dimension",
    "unit": None,
    "description": "Number of data points in the Y-dimension (often the indirectly acquired dimension or F1)."
  },
  "RES": { # Often same as SSX for Bruker
    "fair_name": "number_of_points_x_dimension_alias_res",
    "unit": None,
    "description": "Number of points in the X-dimension (Bruker 'RES' parameter, typically equals SSX)."
  },
  "REY": { # Often same as SSY for Bruker
    "fair_name": "number_of_points_y_dimension_alias_rey",
    "unit": None,
    "description": "Number of points in the Y-dimension (Bruker 'REY' parameter, typically equals SSY)."
  },

  # --- X-Dimension (often Magnetic Field) Parameters ---
  "JEX": {
    "fair_name": "x_dimension_experiment_type_jcamp",
    "unit": None,
    "description": "JCAMP-DX experiment type for the X-dimension (e.g., 'field-sweep')."
  },
  "XXLB": {
    "fair_name": "x_dimension_axis_start_value",
    "unit": None, # Unit defined by XXUN
    "description": "Start value for the X-dimension axis."
  },
  "XXWI": {
    "fair_name": "x_dimension_axis_sweep_width",
    "unit": None, # Unit defined by XXUN
    "description": "Sweep width for the X-dimension axis."
  },
  "XXUN": {
    "fair_name": "x_dimension_axis_unit",
    "unit": None, # The value itself is the unit
    "description": "Unit for the X-dimension axis (e.g., 'G' for Gauss)."
  },
  "HCF": {
    "fair_name": "x_dimension_magnetic_field_center",
    "unit": "G",
    "description": "Center of the magnetic field sweep (typically for the X-dimension)."
  },
  "HSW": {
    "fair_name": "x_dimension_magnetic_field_sweep_width",
    "unit": "G",
    "description": "Width of the magnetic field sweep (typically for the X-dimension)."
  },
  "JUN": { # Primarily for 1D, but might appear if X-dim is the main one
    "fair_name": "primary_x_axis_unit_jcamp",
    "unit": None, # The value itself is the unit, e.g., 'G'
    "description": "JCAMP-DX unit for the primary X-axis (often the X-dimension in 2D, e.g., 'G' for Gauss)."
  },

  # --- Y-Dimension (e.g., MW Power, Time, etc.) Parameters ---
  "JEY": {
    "fair_name": "y_dimension_experiment_type_jcamp",
    "unit": None,
    "description": "JCAMP-DX experiment type for the Y-dimension (e.g., 'mw-power-sweep')."
  },
  "XYLB": {
    "fair_name": "y_dimension_axis_start_value",
    "unit": None, # Unit defined by XYUN
    "description": "Start value for the Y-dimension axis."
  },
  "XYWI": {
    "fair_name": "y_dimension_axis_sweep_width",
    "unit": None, # Unit defined by XYUN
    "description": "Sweep width for the Y-dimension axis."
  },
  "XYUN": {
    "fair_name": "y_dimension_axis_unit",
    "unit": None, # The value itself is the unit
    "description": "Unit for the Y-dimension axis (e.g., 'dB' for decibels)."
  },
  "MPS": {
    "fair_name": "y_dimension_sweep_step_size",
    "unit": None, # Unit typically same as XYUN or implied by context (e.g., dB for power)
    "description": "Step size or increment for the Y-dimension sweep (e.g., microwave power step in dB)."
  },

  # --- Spectrometer / Acquisition Parameters (often apply to X-dim or per point) ---
  "GST": { # Display g-factor for X-axis if field sweep
    "fair_name": "g_factor_start_display",
    "unit": None, # Dimensionless
    "description": "Start value for g-factor axis display/calculation (if X-axis is magnetic field)."
  },
  "GSI": { # Display g-factor for X-axis if field sweep
    "fair_name": "g_factor_increment_display",
    "unit": None, # Dimensionless
    "description": "Increment value for g-factor axis display/calculation (if X-axis is magnetic field)."
  },
  "JNS": {
    "fair_name": "number_of_scans",
    "unit": None, # Integer count
    "description": "Number of accumulated scans (for each point in Y-dimension in a 2D experiment, or total for 1D)."
  },
  "EMF": {
    "fair_name": "magnetic_field_modulation_frequency",
    "unit": "Hz",
    "description": "Frequency of the magnetic field modulation."
  },
  "RCT": {
    "fair_name": "receiver_conversion_time",
    "unit": "ms",
    "description": "Analog-to-digital conversion time per data point."
  },
  "RTC": {
    "fair_name": "receiver_time_constant",
    "unit": "ms",
    "description": "Time constant of the receiver signal filter."
  },
  "RRG": {
    "fair_name": "receiver_gain",
    "unit": None, # Often a relative value or code, sometimes dB
    "description": "Overall gain of the receiver (often a relative value or code)."
  },
  "ROF": {
    "fair_name": "receiver_offset",
    "unit": "%", # Or could be ADC units depending on system
    "description": "Receiver DC offset, often as a percentage of full scale or in ADC units."
  },
  "RMA": {
    "fair_name": "magnetic_field_modulation_amplitude",
    "unit": "G", # Typically peak-to-peak
    "description": "Amplitude of the magnetic field modulation (typically peak-to-peak)."
  },
  "RPH": {
    "fair_name": "receiver_detection_phase",
    "unit": "°",
    "description": "Phase of the phase-sensitive detector (PSD)."
  },
  "RHA": {
    "fair_name": "receiver_detection_harmonic",
    "unit": None, # Integer, e.g., 1 for fundamental, 2 for second harmonic
    "description": "Harmonic of the modulation frequency used for detection (e.g., 1 or 2)."
  },
  "RMF": {
    "fair_name": "receiver_microwave_attenuation", # Or "Reference Arm Microwave Factor"
    "unit": "dB",
    "description": "Microwave attenuation in the receiver path or reference arm (ESP specific)."
  },
  "MF": {
    "fair_name": "microwave_frequency",
    "unit": "GHz",
    "description": "Microwave frequency (ESP)."
  },
  "MPD": {
    "fair_name": "microwave_power_attenuation_setting_initial",
    "unit": "dB",
    "description": "Microwave power attenuation setting in dB, relative to maximum source power. For power sweeps (2D), this is often the initial attenuation value for the Y-dimension."
  }
}

# --- Helper Functions for Saving ---

def _process_parameters(pars: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Processes the raw parameters dictionary using the BRUKER_PARAM_MAP."""
    fair_metadata = {}
    unmapped_parameters = {}

    # Add explicit conversion metadata
    fair_metadata["conversion_info"] = {
        "value": { # Nest value to match other parameter structure
            "converter_script": "fair_converter.py", # Or add versioning
            "conversion_timestamp": datetime.now().isoformat(),
            "eprload_version": "unknown" # TODO: Add version to eprload if possible
        },
        "unit": "",
        "description": "Information about the conversion process to FAIR format."
    }


    for key, value in pars.items():
        # Skip internal keys added by eprload (like abscissa units)
        if key.startswith('_'):
            continue

        if key in BRUKER_PARAM_MAP:
            map_info = BRUKER_PARAM_MAP[key]
            fair_key = map_info["fair_name"]

            # Determine unit: Use map unit, but check if specific axis unit overrides it
            unit = map_info["unit"]
            if unit == "determined_by_XUNI" and "XUNI" in pars:
                unit = pars.get("XUNI", "unknown")
            elif unit == "determined_by_YUNI" and "YUNI" in pars:
                 # Handle potential list format for YUNI/ZUNI if IKKF indicates multiple traces
                 yuni_val = pars.get("YUNI", "unknown")
                 unit = yuni_val.split(',')[0].strip() if isinstance(yuni_val, str) else str(yuni_val)

            elif unit == "determined_by_ZUNI" and "ZUNI" in pars:
                 zuni_val = pars.get("ZUNI", "unknown")
                 unit = zuni_val.split(',')[0].strip() if isinstance(zuni_val, str) else str(zuni_val)


            fair_metadata[fair_key] = {
                "value": value,
                "unit": unit,
                "description": map_info["description"]
            }
        else:
            # Store unmapped parameters separately
            unmapped_parameters[key] = value
            print(f"Parameter '{key}' not found in BRUKER_PARAM_MAP. Storing in 'unmapped_parameters'.")

    return fair_metadata, unmapped_parameters


def _save_to_csv_json(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str
) -> None:
    """Saves data to CSV and structured metadata to JSON."""
    json_file = output_basename.with_suffix(".json")
    csv_file = output_basename.with_suffix(".csv")

    fair_meta, unmapped_meta = _process_parameters(pars)

    print(f"  Saving structured metadata to: {json_file}")
    # --- Save Metadata to JSON ---
    metadata_to_save = {
        "original_file": original_file_path,
        "fair_metadata": fair_meta,
    }
    if unmapped_meta: # Only add if there are unmapped parameters
        metadata_to_save["unmapped_parameters"] = unmapped_meta

    try:
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(metadata_to_save, f, indent=4, default=str) # default=str handles non-serializables
    except IOError as e:
        warnings.warn(f"Could not write JSON file {json_file}: {e}")
    except TypeError as e:
        warnings.warn(f"Error serializing metadata to JSON for {json_file}: {e}. Some parameters might not be saved correctly.")

    print(f"  Saving data to: {csv_file}")
    # --- Save Data to CSV ---
    # (CSV saving logic remains the same as before, it doesn't use the detailed map directly)
    try:
        with open(csv_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)

            # --- Write Header ---
            writer.writerow(['# EPR Data Export'])
            writer.writerow(['# Original File:', original_file_path])
            # Add a few key parameters derived from the FAIR map if available
            mwfq_info = fair_meta.get('microwave_frequency', {})
            field_info = fair_meta.get('field_center', fair_meta.get('field_sweep_start', {}))
            sweep_info = fair_meta.get('field_sweep_width', fair_meta.get('field_sweep_increment', {}))

            writer.writerow(['# Microwave_Frequency:', f"{mwfq_info.get('value', 'N/A')} {mwfq_info.get('unit', '')}".strip()])
            writer.writerow(['# Field_Center/Start:', f"{field_info.get('value', 'N/A')} {field_info.get('unit', '')}".strip()])
            writer.writerow(['# Field_Sweep/Increment:', f"{sweep_info.get('value', 'N/A')} {sweep_info.get('unit', '')}".strip()])
            writer.writerow(['# Data_Shape:', str(y.shape)])
            writer.writerow(['# Data_Type:', str(y.dtype)])
            writer.writerow(['# ---'])

            # --- Prepare Data Columns ---
            header_row = []
            data_columns = []
            is_complex = np.iscomplexobj(y)
            is_2d = y.ndim == 2

            # Get Abscissa Units using potentially mapped keys first
            x_unit_val = fair_meta.get('axis_x_unit', {}).get('value', 'a.u.')
            # Handle list case from IKKF in BES3T for Y/Z units if needed
            y_unit_val = fair_meta.get('axis_y_unit', {}).get('value', 'a.u.')
            if isinstance(y_unit_val, str) and ',' in y_unit_val: y_unit_val = y_unit_val.split(',')[0].strip()
            z_unit_val = fair_meta.get('axis_z_unit', {}).get('value', 'a.u.')
            if isinstance(z_unit_val, str) and ',' in z_unit_val: z_unit_val = z_unit_val.split(',')[0].strip()


            if not is_2d: # 1D Data
                n_pts = y.shape[0]
                # Abscissa Column
                if x is not None and isinstance(x, np.ndarray) and x.shape == y.shape:
                    header_row.append(f"Abscissa ({x_unit_val})")
                    data_columns.append(x)
                else:
                    header_row.append("Index")
                    data_columns.append(np.arange(n_pts))
                    if x is not None: warnings.warn("Provided x-axis ignored for CSV (shape mismatch or not ndarray). Using index.")

                # Intensity Column(s)
                if is_complex:
                    header_row.extend(["Intensity_Real (a.u.)", "Intensity_Imag (a.u.)"])
                    data_columns.append(np.real(y))
                    data_columns.append(np.imag(y))
                else:
                    header_row.append("Intensity (a.u.)")
                    data_columns.append(y)

            else: # 2D Data - use "long" format (X, Y, Value(s))
                ny, nx = y.shape
                x_coords_flat = np.arange(nx) # Default: Index
                y_coords_flat = np.arange(ny) # Default: Index
                header_row.extend([f"X_Index ({nx} points)", f"Y_Index ({ny} points)"])

                # Determine X and Y axes from input 'x'
                if isinstance(x, list) and len(x) >= 2:
                    x_axis, y_axis = x[0], x[1]
                    if isinstance(x_axis, np.ndarray) and x_axis.size == nx:
                        x_coords_flat = x_axis
                        header_row[0] = f"X_Axis ({x_unit_val})"
                    if isinstance(y_axis, np.ndarray) and y_axis.size == ny:
                        y_coords_flat = y_axis
                        header_row[1] = f"Y_Axis ({y_unit_val})" # Use Y unit here
                elif isinstance(x, np.ndarray) and x.ndim == 1 and x.size == nx:
                     x_coords_flat = x
                     header_row[0] = f"X_Axis ({x_unit_val})"

                # Create grid and flatten
                xx, yy = np.meshgrid(x_coords_flat, y_coords_flat)
                data_columns.append(xx.ravel())
                data_columns.append(yy.ravel())

                # Intensity Column(s)
                if is_complex:
                    header_row.extend(["Intensity_Real (a.u.)", "Intensity_Imag (a.u.)"])
                    data_columns.append(np.real(y).ravel())
                    data_columns.append(np.imag(y).ravel())
                else:
                    header_row.append("Intensity (a.u.)")
                    data_columns.append(y.ravel())

            # --- Write Data ---
            writer.writerow(header_row)
            rows_to_write = np.stack(data_columns, axis=-1)
            writer.writerows(rows_to_write)

    except IOError as e:
        warnings.warn(f"Could not write CSV file {csv_file}: {e}")
    except Exception as e:
        warnings.warn(f"An unexpected error occurred while writing CSV {csv_file}: {e}")

def _try_set_h5_attr(h5_object, key: str, value: Any):
    """Helper to safely set HDF5 attributes, converting to string on type error."""
    try:
        if value is None:
            h5_object.attrs[key] = "None" # Store None as string
        elif isinstance(value, (list, tuple)) and all(isinstance(i, (int, float, str, np.number, bytes)) for i in value):
             # Store lists/tuples of basic types directly if h5py supports it, else convert
             try:
                  h5_object.attrs[key] = value # Let h5py handle if possible (newer versions are better)
             except TypeError:
                 # Try converting to numpy array first
                 try:
                     h5_object.attrs[key] = np.array(value)
                 except TypeError:
                      h5_object.attrs[key] = str(value) # Fallback to string
        elif isinstance(value, Path): # Convert Path object to string
            h5_object.attrs[key] = str(value)
        else:
            h5_object.attrs[key] = value # Try direct assignment

    except TypeError:
        # If type is not directly supported by HDF5 attributes, convert to string
        warnings.warn(f"Could not store attribute '{key}' (type: {type(value)}) directly in HDF5 attributes. Converting to string.")
        h5_object.attrs[key] = str(value)
    except Exception as e:
         warnings.warn(f"Unexpected error storing attribute '{key}': {type(e).__name__} - {e}. Skipping.")


def _save_to_hdf5(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str
) -> None:
    """Saves data and structured metadata to an HDF5 file."""
    h5_file = output_basename.with_suffix(".h5")
    fair_meta, unmapped_meta = _process_parameters(pars)

    print(f"  Saving structured data and metadata to: {h5_file}")

    try:
        with h5py.File(h5_file, 'w') as f:
            # --- Store Global Metadata ---
            f.attrs["original_file"] = original_file_path
            f.attrs["description"] = "FAIR representation of EPR data converted from Bruker format."
            f.attrs["conversion_timestamp"] = datetime.now().isoformat()
            f.attrs["converter_script_version"] = "fair_converter_v1.1" # Example version

            # --- Store Structured FAIR Metadata ---
            param_grp = f.create_group("metadata/parameters_fair")
            param_grp.attrs["description"] = "Mapped parameters with units and descriptions."

            for fair_key, info in fair_meta.items():
                item_grp = param_grp.create_group(fair_key)
                _try_set_h5_attr(item_grp, "value", info["value"])
                _try_set_h5_attr(item_grp, "unit", info["unit"])
                _try_set_h5_attr(item_grp, "description", info["description"])

            # --- Store Unmapped Parameters ---
            if unmapped_meta:
                unmap_grp = f.create_group("metadata/parameters_original")
                unmap_grp.attrs["description"] = "Parameters from the original file not found in the FAIR mapping."
                for key, value in unmapped_meta.items():
                    _try_set_h5_attr(unmap_grp, key, value)

            # --- Store Data ---
            data_grp = f.create_group("data")
            ds_y = data_grp.create_dataset("intensity", data=y)
            ds_y.attrs["description"] = "Experimental intensity data."
            ds_y.attrs["units"] = "a.u."
            if np.iscomplexobj(y):
                 ds_y.attrs["signal_type"] = "complex"
            else:
                 ds_y.attrs["signal_type"] = "real"

            # Get Abscissa Units from FAIR map
            x_unit_val = fair_meta.get('axis_x_unit', {}).get('value', 'a.u.')
            y_unit_val = fair_meta.get('axis_y_unit', {}).get('value', 'a.u.')
            if isinstance(y_unit_val, str) and ',' in y_unit_val: y_unit_val = y_unit_val.split(',')[0].strip()
            z_unit_val = fair_meta.get('axis_z_unit', {}).get('value', 'a.u.')
            if isinstance(z_unit_val, str) and ',' in z_unit_val: z_unit_val = z_unit_val.split(',')[0].strip()

            # Store abscissa(e)
            axis_datasets = {} # Store references to axis datasets
            if x is None:
                if y.ndim >= 1:
                    nx = y.shape[-1]
                    ds_x = data_grp.create_dataset("abscissa_x", data=np.arange(nx))
                    _try_set_h5_attr(ds_x, "units", "points")
                    _try_set_h5_attr(ds_x, "description", "X axis (index)")
                    _try_set_h5_attr(ds_x, "axis_type", "index")
                    axis_datasets['x'] = ds_x
                if y.ndim >= 2:
                    ny = y.shape[-2]
                    ds_y_ax = data_grp.create_dataset("abscissa_y", data=np.arange(ny))
                    _try_set_h5_attr(ds_y_ax, "units", "points")
                    _try_set_h5_attr(ds_y_ax, "description", "Y axis (index)")
                    _try_set_h5_attr(ds_y_ax, "axis_type", "index")
                    axis_datasets['y'] = ds_y_ax
                # Add Z if needed

            elif isinstance(x, np.ndarray): # 1D data
                ds_x = data_grp.create_dataset("abscissa_x", data=x)
                _try_set_h5_attr(ds_x, "units", x_unit_val)
                _try_set_h5_attr(ds_x, "description", f"X axis ({fair_meta.get('axis_x_unit',{}).get('description','Unknown')})")
                _try_set_h5_attr(ds_x, "axis_type", "independent_variable")
                axis_datasets['x'] = ds_x
            elif isinstance(x, list): # Multi-D data
                if len(x) >= 1 and x[0] is not None and isinstance(x[0], np.ndarray):
                    ds_x = data_grp.create_dataset("abscissa_x", data=x[0])
                    _try_set_h5_attr(ds_x, "units", x_unit_val)
                    _try_set_h5_attr(ds_x, "description", f"X axis ({fair_meta.get('axis_x_unit',{}).get('description','Unknown')})")
                    _try_set_h5_attr(ds_x, "axis_type", "independent_variable_x")
                    axis_datasets['x'] = ds_x
                if len(x) >= 2 and x[1] is not None and isinstance(x[1], np.ndarray):
                     ds_y_ax = data_grp.create_dataset("abscissa_y", data=x[1])
                     _try_set_h5_attr(ds_y_ax, "units", y_unit_val)
                     _try_set_h5_attr(ds_y_ax, "description", f"Y axis ({fair_meta.get('axis_y_unit',{}).get('description','Unknown')})")
                     _try_set_h5_attr(ds_y_ax, "axis_type", "independent_variable_y")
                     axis_datasets['y'] = ds_y_ax
                if len(x) >= 3 and x[2] is not None and isinstance(x[2], np.ndarray):
                     ds_z_ax = data_grp.create_dataset("abscissa_z", data=x[2])
                     _try_set_h5_attr(ds_z_ax, "units", z_unit_val)
                     _try_set_h5_attr(ds_z_ax, "description", f"Z axis ({fair_meta.get('axis_z_unit',{}).get('description','Unknown')})")
                     _try_set_h5_attr(ds_z_ax, "axis_type", "independent_variable_z")
                     axis_datasets['z'] = ds_z_ax

            # --- Link axes to data dimensions using HDF5 Dimension Scales API ---
            # Use explicit positive indexing and add error handling
            if 'intensity' in data_grp:
                dims = ds_y.dims
                current_ndim = ds_y.ndim # Use ndim from the actual HDF5 dataset

                # Link X dimension (last dimension, index ndim-1)
                if current_ndim >= 1 and 'x' in axis_datasets:
                    x_dim_index = current_ndim - 1
                    try:
                        dims[x_dim_index].label = 'x'
                        dims[x_dim_index].attach_scale(axis_datasets['x'])
                    except OverflowError as e_ovfl:
                         warnings.warn(f"OverflowError setting label/scale for X dim (index {x_dim_index}): {e_ovfl}. Skipping dimension link.")
                    except Exception as e:
                         warnings.warn(f"Error setting label/scale for X dim (index {x_dim_index}): {type(e).__name__} - {e}. Skipping dimension link.")

                # Link Y dimension (second to last dimension, index ndim-2)
                if current_ndim >= 2 and 'y' in axis_datasets:
                    y_dim_index = current_ndim - 2
                    try:
                        dims[y_dim_index].label = 'y'
                        dims[y_dim_index].attach_scale(axis_datasets['y'])
                    except OverflowError as e_ovfl:
                         warnings.warn(f"OverflowError setting label/scale for Y dim (index {y_dim_index}): {e_ovfl}. Skipping dimension link.")
                    except Exception as e:
                         warnings.warn(f"Error setting label/scale for Y dim (index {y_dim_index}): {type(e).__name__} - {e}. Skipping dimension link.")

                # Link Z dimension (third to last dimension, index ndim-3)
                if current_ndim >= 3 and 'z' in axis_datasets:
                    z_dim_index = current_ndim - 3
                    try:
                        dims[z_dim_index].label = 'z'
                        dims[z_dim_index].attach_scale(axis_datasets['z'])
                    except OverflowError as e_ovfl:
                         warnings.warn(f"OverflowError setting label/scale for Z dim (index {z_dim_index}): {e_ovfl}. Skipping dimension link.")
                    except Exception as e:
                         warnings.warn(f"Error setting label/scale for Z dim (index {z_dim_index}): {type(e).__name__} - {e}. Skipping dimension link.")


    except ImportError:
        warnings.warn("h5py library not found. Skipping HDF5 output. Install with 'pip install h5py'")
    except IOError as e:
        warnings.warn(f"Could not write HDF5 file {h5_file}: {e}")
    except Exception as e:
         warnings.warn(f"An unexpected error occurred while writing HDF5 file {h5_file}: {type(e).__name__} - {e}")
         import traceback
         traceback.print_exc()

# --- Main Conversion Function ---

def append_fair_metadata(pars: Dict[str, Any]) -> Dict[str, Any]:
    """
    Appends FAIR metadata (using BRUKER_PARAM_MAP) to the input parameter dictionary `pars`
    under the key '_fair_metadata'. Returns the updated dictionary.

    Args:
        pars (dict): Original parameter dictionary.

    Returns:
        dict: Updated parameter dictionary with '_fair_metadata' key added.
    """
    fair_metadata, unmapped = _process_parameters(pars)
    pars['_fair_metadata'] = fair_metadata
    pars['_unmapped_parameters'] = unmapped
    return pars

def convert_bruker_to_fair(
    input_file_or_dir: Union[str, Path, None] = None,
    output_dir: Optional[Union[str, Path]] = None,
    scaling: str = "",
    output_formats: List[str] = ['csv_json', 'hdf5']
) -> None:
    """
    Loads Bruker EPR data using eprload and converts it to FAIR formats
    with structured metadata.

    Args:
        input_file_or_dir (Union[str, Path, None], optional): Path to the Bruker
            data file (.dta, .dsc, .spc, .par) or a directory. If None or
            a directory, a file dialog will open. Defaults to None.
        output_dir (Optional[Union[str, Path]], optional): Directory to save
            the converted files. If None, files are saved in the same
            directory as the input file. Defaults to None.
        scaling (str, optional): Scaling options passed to eprload
            (e.g., 'nPGT'). Defaults to "".
        output_formats (List[str], optional): List of formats to generate.
            Options: 'csv_json', 'hdf5'. Defaults to ['csv_json', 'hdf5'].

    Returns:
        None. Prints status messages and warnings.
    """
    print(f"Starting FAIR conversion process...")

    # Load data using eprload (suppress its internal plotting)
    x, y, pars, original_file_path_str = eprload(
        file_name=input_file_or_dir,
        scaling=scaling,
        plot_if_possible=False # Important: disable plotting in eprload
    )

    # --- Check if loading was successful ---
    if y is None or pars is None or original_file_path_str is None:
        print("Data loading failed or was cancelled by user. Aborting conversion.")
        return

    print(f"Successfully loaded: {original_file_path_str}")
    original_file_path = Path(original_file_path_str)

    # --- Determine Output Location and Basename ---
    if output_dir is None:
        output_path = original_file_path.parent
    else:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True) # Create output dir if needed

    output_basename = output_path / original_file_path.stem
    print(f"Output base name: {output_basename}")

    # --- Perform Conversions based on requested formats ---
    if not output_formats:
         warnings.warn("No output formats specified. Nothing to do.")
         return

    print("\nProcessing parameters and generating outputs...")
    # Note: _process_parameters is called inside each save function now

    if 'csv_json' in output_formats:
        _save_to_csv_json(output_basename, x, y, pars, original_file_path_str)
        print("...CSV + JSON generation finished.")

    if 'hdf5' in output_formats:
         _save_to_hdf5(output_basename, x, y, pars, original_file_path_str)
         print("...HDF5 generation finished.")

    print("\nFAIR conversion process finished.")


# --- Example Usage ---
def save_fair(
    output_basename: Union[str, Path],
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Dict[str, Any],
    original_file_path: str,
    output_formats: List[str] = ['csv_json', 'hdf5']
) -> None:
    """
    Saves already-loaded EPR data to one or more FAIR formats.

    Args:
        output_basename (Union[str, Path]): The base path for the output files,
            without an extension.
        x (Union[np.ndarray, List[np.ndarray], None]): Abscissa data.
        y (np.ndarray): Ordinate data.
        params (Dict[str, Any]): Dictionary of parameters.
        original_file_path (str): The full path of the original loaded file.
        output_formats (List[str], optional): A list of formats to save.
            Defaults to ['csv_json', 'hdf5'].
    """
    output_basename = Path(output_basename)
    output_path = output_basename.parent
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Saving data from '{original_file_path}' to FAIR formats...")

    if 'csv_json' in output_formats:
        _save_to_csv_json(output_basename, x, y, params, original_file_path)
        print("...CSV + JSON generation finished.")

    if 'hdf5' in output_formats:
        _save_to_hdf5(output_basename, x, y, params, original_file_path)
        print("...HDF5 generation finished.")

    print("\nFAIR saving process finished.")


if __name__ == "__main__":
    print("--- Running Bruker to FAIR Converter Example (with structured metadata) ---")
    print("This will first use eprload's dialog to select a Bruker file,")
    print("then convert it to CSV/JSON and HDF5 formats in the same directory,")
    print("applying the parameter map for richer metadata.")
    print("-" * 50)

    # Example 1: Use file dialog to select input, save to same directory
    convert_bruker_to_fair()

    # # Example 2: Specify an input file and output directory
    # input_bruker_file = "path/to/your/test_data.DSC" # <-- CHANGE THIS
    # output_folder = "path/to/your/output_folder"   # <-- CHANGE THIS
    # input_path = Path(input_bruker_file)
    #
    # if input_path.exists():
    #     print("\n" + "-" * 50)
    #     print(f"--- Example 2: Processing specified file: {input_bruker_file} ---")
    #     convert_bruker_to_fair(
    #         input_file_or_dir=input_bruker_file,
    #         output_dir=output_folder,
    #         scaling="nG" # Example scaling
    #     )
    # else:
    #      print(f"\nSkipping Example 2: Input file not found at '{input_bruker_file}'")

    print("-" * 50)
    print("--- Converter Example Finished ---")