#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1024
XMIN	100.000000
XWID	6000.000000
*
*	Documentational Text:
*
TITL	'Er'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    04/13/06
TIME    09:46:01
CMNT
SAMP    NTE 36
SFOR    bois pyrolyse
STAG    A
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3
A1CT    0.31
A1SW    0.6
MWFQ    9.704197e+09
MWPW    1.013e-05
AVGS    0
RESO    md4w10301
SPTP    0.08192
RCAG    60
RCHM    1
B0MA    0.0001
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    1024
RCTC    0.04096
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     fieldCtrl, 1.0

CenterField        3100.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
SweepDirection     Up
SweepWidth         6000.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.704197 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
Power              0.01013 mW
PowerAtten         43 dB

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         0
NbScansDone        0
NbScansToDo        1
ReplaceMode        Off

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
Calibrated         True
ConvTime           81.92 ms
DModAFCTrap        True
DModAmp            1.00 G
DModCalibrated     True
DModDetectSCT      First
DModEliDelay       1.0 us
DModExtLockIn      False
DModExtTrigger     False
DModFieldMod       First
DModGain           60 dB
DModHighPass       True
DModIntegrator     True
DModModOutput      Internal
DModSignalInput    Internal
DModTimeConst      1.28 ms
DoubleModFreq      5.00 kHz
DoubleModPhase     0.0
DoubleMode         False
EliDelay           1.0 us
ExtLockIn          False
ExtTrigger         False
Gain               60 dB
Harmonic           1
HighPass           True
Integrator         True
ModAmp             1.00 G
ModFreq            100.00 kHz
ModInput           Internal
ModOutput          Internal
ModPhase           0.0
Offset             0.0 %
QuadMode           False
QuadPhase          90.0
Resolution         1024
Resonator          1
SctNorm            False
SctRevision        SCT
SignalInput        Internal
SweepTime          83.89 s
TimeConst          40.96 ms
TuneCaps           32

*
************************************************************
