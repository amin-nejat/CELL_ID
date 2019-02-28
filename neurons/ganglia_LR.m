% Create a list of ganglia & their neurons.
ganglion_LR = {};

% Anterior Pharyngeal Bulb
ganglion_LR(1).name = 'Anterior Pharyngeal Bulb';
ganglion_LR(1).neurons = {
    'I1L'
    'I1R'
    'I2L'
    'I2R'
    'I3'
    'M3L'
    'M3R'
    'M4'
    'MCL'
    'MCR'
    'MI'
    'NSML'
    'NSMR'
    };

% Posterior Pharyngeal Bulb
ganglion_LR(end+1).name = 'Posterior Pharyngeal Bulb';
ganglion_LR(end).neurons = {
    'I4'
    'I5'
    'I6'
    'M1'
    'M2L'
    'M2R'
    'M5'
    };

% Anterior Ganglion (Left)
ganglion_LR(end+1).name = 'Anterior Ganglion (Left)';
ganglion_LR(end).neurons = {
    'BAGL'
    'CEPVL'
    'IL1L'
    'IL1DL'
    'IL1VL'
    'IL2L'
    'IL2DL'
    'IL2VL'
    'OLLL'
    'OLQDL'
    'OLQVL'
    'RIPL'
    'RMEL'
    'RMED'
    'RMEV'
    'URADL'
    'URAVL'
    'URBL'
    'URYDL'
    'URYVL'
    };

% Anterior Ganglion (Right)
ganglion_LR(end+1).name = 'Anterior Ganglion (Right)';
ganglion_LR(end).neurons = {
    'BAGR'
    'CEPVR'
    'IL1R'
    'IL1DR'
    'IL1VR'
    'IL2R'
    'IL2DR'
    'IL2VR'
    'OLLR'
    'OLQDR'
    'OLQVR'
    'RIPR'
    'RMER'
    'RMED'
    'RMEV'
    'URADR'
    'URAVR'
    'URBR'
    'URYDR'
    'URYVR'
    };

% Dorsal Ganglion
ganglion_LR(end+1).name = 'Dorsal Ganglion';
ganglion_LR(end).neurons = {
    'ALA'
    'CEPDL'
    'CEPDR'
    'RID'
    'URXL'
    'URXR'
    };

% Lateral Ganglion (Left)
ganglion_LR(end+1).name = 'Lateral Ganglion (Left)';
ganglion_LR(end).neurons = {
    'ADAL'
    'ADEL'
    'ADFL'
    'ADLL'
    'AFDL'
    'AIBL'
    'AINL'
    'AIZL'
    'ASEL'
    'ASGL'
    'ASHL'
    'ASIL'
    'ASJL'
    'ASKL'
    'AUAL'
    'AVAL'
    'AVBL'
    'AVDL'
    'AVEL'
    'AVHL'
    'AVJL'
    'AWAL'
    'AWBL'
    'AWCL'
    'AWCL(OFF)'
    'AWCL(ON)'
    'FLPL'
    'RIAL'
    'RIBL'
    'RICL'
    'RIML'
    'RIVL'
    'RMDL'
    'RMDVL'
    'RMGL'
    'SAAVL'
    'SIBDL'
    'SMDVL'
    };

% Lateral Ganglion (Right)
ganglion_LR(end+1).name = 'Lateral Ganglion (Right)';
ganglion_LR(end).neurons = {
    'ADAR'
    'ADER'
    'ADFR'
    'ADLR'
    'AFDR'
    'AIBR'
    'AINR'
    'AIZR'
    'AQR'
    'ASER'
    'ASGR'
    'ASHR'
    'ASIR'
    'ASJR'
    'ASKR'
    'AUAR'
    'AVAR'
    'AVBR'
    'AVDR'
    'AVER'
    'AVHR'
    'AVJR'
    'AWAR'
    'AWBR'
    'AWCR'
    'AWCR(OFF)'
    'AWCR(ON)'
    'FLPR'
    'RIAR'
    'RIBR'
    'RICR'
    'RIMR'
    'RIVR'
    'RMDR'
    'RMDVR'
    'RMGR'
    'SAAVR'
    'SIBDR'
    'SMDVR'
    };

% Ventral Ganglion
ganglion_LR(end+1).name = 'Ventral Ganglion';
ganglion_LR(end).neurons = {
    'AIAL'
    'AIAR'
    'AIML'
    'AIMR'
    'AIYL'
    'AIYR'
    'AVKL'
    'AVKR'
    'AVL'
    'SAADL'
    'SAADR'
    'SIADL'
    'SIADR'
    'SIAVL'
    'SIAVR'
    'SIBVL'
    'SIBVR'
    'SMBDL'
    'SMBDR'
    'SMBVL'
    'SMBVR'
    'SMDDL'
    'SMDDR'
    'RIH'
    'RIR'
    'RIS'
    'RMDDL'
    'RMDDR'
    'RMFL'
    'RMFR'
    'RMHL'
    'RMHR'
    };

% Retrovesicular Ganglion
ganglion_LR(end+1).name = 'Retrovesicular Ganglion';
ganglion_LR(end).neurons = {
    'AS01'
    'AVFL'
    'AVFR'
    'AVG'
    'DA01'
    'DB01'
    'DB02'
    'DD01'
    'RIFL'
    'RIFR'
    'RIGL'
    'RIGR'
    'SABD'
    'SABVL'
    'SABVR'
    'VA01'
    'VD01'
    'VD02'
    'VB01'
    'VB02'
    };

% Midbody (Left)
ganglion_LR(end+1).name = 'Midbody (Left)';
ganglion_LR(end).neurons = {
    'ALML'
    'BDUL'
    'CANL'
    'HSNL'
    'PDEL'
    'PVDL'
    'PVM'
    'SDQL'
    };

% Midbody (Right)
ganglion_LR(end+1).name = 'Midbody (Right)';
ganglion_LR(end).neurons = {
    'ALMR'
    'AVM'
    'BDUR'
    'CANR'
    'HSNR'
    'PDER'
    'PVDR'
    'SDQR'
    };

% Ventral Nerve Cord
ganglion_LR(end+1).name = 'Ventral Nerve Cord';
ganglion_LR(end).neurons = {
    'AS02'
    'AS03'
    'AS04'
    'AS05'
    'AS06'
    'AS07'
    'AS08'
    'AS09'
    'AS10'
    'DA02'
    'DA03'
    'DA04'
    'DA05'
    'DA06'
    'DA07'
    'DB03'
    'DB04'
    'DB05'
    'DB06'
    'DB07'
    'DD02'
    'DD03'
    'DD04'
    'DD05'
    'VA02'
    'VA03'
    'VA04'
    'VA05'
    'VA06'
    'VA07'
    'VA08'
    'VA09'
    'VA10'
    'VA11'
    'VB03'
    'VB04'
    'VB05'
    'VB06'
    'VB07'
    'VB08'
    'VB09'
    'VB10'
    'VB11'
    'VC01'
    'VC02'
    'VC03'
    'VC04'
    'VC05'
    'VC06'
    'VD03'
    'VD04'
    'VD05'
    'VD06'
    'VD07'
    'VD08'
    'VD09'
    'VD10'
    'VD11'
    };

% Pre-Anal Ganglion
ganglion_LR(end+1).name = 'Pre-Anal Ganglion';
ganglion_LR(end).neurons = {
    'AS11'
    'DA08'
    'DA09'
    'DD06'
    'PDA'
    'PDB'
    'PVPL'
    'PVPR'
    'PVT'
    'VA12'
    'VD12'
    'VD13'
    };

% Dorso-Rectal Ganglion
ganglion_LR(end+1).name = 'Dorso-Rectal Ganglion';
ganglion_LR(end).neurons = {
    'DVA'
    'DVB'
    'DVC'
    };

% Tail (Left)
ganglion_LR(end+1).name = 'Tail (Left)';
ganglion_LR(end).neurons = {
    'ALNL'
    'LUAL'
    'PHAL'
    'PHBL'
    'PHCL'
    'PLML'
    'PLNL'
    'PQR'
    'PVCL'
    'PVNL'
    'PVQL'
    'PVWL'
    };

% Tail (Right)
ganglion_LR(end+1).name = 'Tail (Right)';
ganglion_LR(end).neurons = {
    'ALNR'
    'LUAR'
    'PHAR'
    'PHBR'
    'PHCR'
    'PLMR'
    'PLNR'
    'PVCR'
    'PVNR'
    'PVQR'
    'PVR'
    'PVWR'
    };

% Save the ganglia.
save('ganglia_LR.mat', 'ganglion_LR');
