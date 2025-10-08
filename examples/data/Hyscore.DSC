#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	IDX
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	256
XMIN	0.000000
XWID	6630.000000
YPTS	256
YMIN	0.000000
YWID	6630.000000
*
*	Documentational Text:
*
TITL	'035_HYSCORE_d1_650_pi_2_36_0deg'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'Time'
YNAM	'Time'
IRUNI	''
IIUNI	''
XUNI	'ns'
YUNI	'ns'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    11/30/17
TIME    04:29:46
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    
OXS1    
AXS1    
AXS2    
AXS3    
MWPW    0.0004743
A1CT    0.348
B0VL    0.3468184
A1SW    0.01
MWFQ    9.727467e+09
AVGS    10
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        38.4
Power              0.4743 mW
PowerAtten         25.0 dB

.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        3480.00 G
Delay              0.0 s
FieldFlyback       On
FieldPosition      3468.184 G
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     3456.221 G
SweepDirection     Up
SweepWidth         100.0 G
WidthTM            200.0 G

.DVC     freqCounter, 1.0

FrequencyMon       9.727467 GHz
QMonitBridge       On

.DVC     ftBridge, 1.0

Attenuation        11.00 dB
ELDORAtt           30 dB
FrequencyA         9.500000 GHz
VideoBW            200 MHz
VideoGain          51 dB

.DVC     ftEpr, 1.0

AWGPhaseShift      90.00 deg.
AutoTimeOut        True
AveragesPerScan    10
FTEzAWGELDORa      100.00 %
FTEzAWGELDORf      0.000 MHz
FTEzAWGELDORw      0.000 MHz
FTEzB0SettlDel     0 s
FTEzB0Step         0 G
FTEzCurrELDOR      9.100000 GHz
FTEzD11deer        0 ns
FTEzD20deer        0 ns
FTEzD3deer         300 ns
FTEzDXdeer         0 ns
FTEzDelay0         350 ns
FTEzDelay1         650 ns
FTEzDelay2         300 ns
FTEzDelay4         0 ns
FTEzDeltaX         26 ns
FTEzDeltaY         26 ns
FTEzELDORStart     9.100000 GHz
FTEzELDORWidth     0.000000 GHz
FTEzENDORDetDel    0 ns
FTEzExtTrgDel      200 ns
FTEzIntgBase       1 ns
FTEzIntgGate       300 ns
FTEzMWPiHalf       36 ns
FTEzNPls           1
FTEzP20deer        0 ns
FTEzP3             40 ns
FTEzRF1Pi          0 ns
FTEzRF1PlsSep      2 ns
FTEzRF1Pump        0.500 MHz
FTEzRF1Start       0.500 MHz
FTEzRF1Width       0.000 MHz
FTEzSRELDORExp     False
FTEzSXdeer         1024
FTEzSXdeerExp      False
FTEzSleep          0 s
FieldIsStatic      False
GradIntPulse       False
GrdEnable          False
LastXAxis          ?
LastYAxis          ?
MDigPrg            {3;400,1,1;1[mixed]} 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
MMWaveLOFreq       26.880000 GHz
OnlyAWGChans       False
PCycleAllowed      True
PCycleOn           True
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     HYSCORE
PlsSPELLISTSlct    4-step
PlsSPELPhPrgEx     Normal
QuadDetect         True
ReplaceMode        Off
ShotRepTime        899.64 us
ShotsPLoop         10
SmoothAllowed      False
SmoothPoints       1
SptProgress        100 %
StochMode          False
SweepsPExp         10
TriggerTimeOut     20 s
XSpecRes           256
YSpecRes           256

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         10
NbScansDone        10
NbScansToDo        10
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       1

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
