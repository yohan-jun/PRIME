function bSuccess = create_pTXRFPulse( afGradient, afRFPulse, varargin)
% %  afRFPulse should be in "V/degree"

% % 	bSuccess = create_pTXRFPulse(dir, filenum, opt())
% % 	
% % 	description
% % 	
% % 	default options:
% %         dopt.RFPULSE_ID         = 0;        % different ID for several pulses used within one sequence 
% %         dopt.DESIGNFLIPANGLE    = 90;
% %         dopt.FACTOROVERSAMPLE   = 1;        % compared to GRADIENT_RASTER_TIME [10 us]
% %         dopt.WRITEOLDFORMAT     = 0;        % option to write old TX format (TX(i)/*.float)
% %         dopt.VERBOSE            = 1;
% % 	
% % 	return: true on Success
% % 	
% % 	Functions called: 
% % 	Tested for: 
% % 	
% % 	June 2009 -  josef.pfeuffer@siemens.com, Kawin Setsompop, MIT
% %
% %     Version history:
% %         'v1.0, 09-06-27 (jp-ks)' first version
% % 	
% % -----------------------------------------------------------------------------
% %   Copyright (C) Siemens AG 2009  All Rights Reserved. Confidential.
% % -----------------------------------------------------------------------------
FCTNAME = 'create_pTXRFPulse'; nargVars = 2;
VERSION = 'v1.0, 09-06-28 (jp-ks)';         % update version info here

PULSE_FILE_NAME = 'pTXRFPulse';             % base output file name

% global defines
bSuccess = 0;           
iFloatDigits = 5;       % accuracy of written floats
          
      %%%%  option/narg handling: default options: define whole struct and comment
dopt.RFPULSE_ID         = 0;        % different ID for several pulses used within one sequence 
dopt.FACTOROVERSAMPLE   = 1;        % compared to GRADIENT_RASTER_TIME [10 us]
dopt.WRITEOLDFORMAT     = 0;        % option to write old TX format (TX(i)/*.float)
dopt.VERBOSE            = 1;
    % RF pulse properties
dopt.PULSENAME          = PULSE_FILE_NAME;
dopt.FAMILY             = 'pTX';
dopt.COMMENT            = 'pTX';
dopt.NOMFLIPANGLE       = 90.;
dopt.AMPLINT            = 100.;
dopt.ABSINT             = 100.;
dopt.POWERINT           = 100.;
dopt.MINSLICE           = 1.;
dopt.MAXSLICE           = 1.;
dopt.REFGRAD            = 1.;

      %%%%  option/narg handling - error: narg > nargVars; demoMode: narg < nargVars
[dopt, narg] 	= handlerOptionsNarg(dopt, varargin, nargin, [0 nargVars]);

      %%%%  option/narg handling: translate options to local parameters
iRFPulseID          = dopt.RFPULSE_ID;
bIsFloatFormat      = dopt.WRITEOLDFORMAT;
iFactorOversample   = dopt.FACTOROVERSAMPLE;
bVerbose            = dopt.VERBOSE;
strPulseName        = dopt.PULSENAME;                   % RF pulse properties 
strFamily           = dopt.FAMILY;          
strComment          = dopt.COMMENT;
fNominalFlipAngle   = dopt.NOMFLIPANGLE;    
fAmplInt            = dopt.AMPLINT;
fAbsInt             = dopt.ABSINT;
fPowerInt           = dopt.POWERINT;
fMinSlice           = dopt.MINSLICE;
fMaxSlice           = dopt.MAXSLICE;
fRefGrad            = dopt.REFGRAD;

% B0B1_DataStructureName 
% RF_DataStructureName  
% PatientB0B1_name
% PatientB0B1_age
% PatientB0B1_weight
% PatientB0B1_position
% PatientB0B1_sex
% SAR_DataStructureName
% SAR_LocalToGlobal
% SAR_CoilType
% SAR_ButlerOrNot
% SAR_Nchn
% TimeStamp
% CheckSum


if (nargin < nargVars); help(FCTNAME); bSuccess = demoMode(FCTNAME); return; end
      %%%%  end: option/narg handling: 

% if requested: write old TX format (TX(i)/*.float)
if( bIsFloatFormat )
    bSuccess = create_header_parallel7T( afGradient, afRFPulse, iFactorOversample );
    if( bVerbose)
        fprintf('%s(): <TX(i)/*.float> successfully written (FloatFormat). <%s>', FCTNAME, VERSION);
    end
    return
end

G_length = length( afGradient );
B_length = G_length*iFactorOversample;
if( rem(length(afRFPulse),B_length) ~= 0 )
    whos afGradient afRFPulse
    error(sprintf('%s: wrong data sizes.', FCTNAME))  
end

iNUsedChannels = length(afRFPulse)/B_length;
for count = 1:iNUsedChannels
    rf(:,count) = afRFPulse(1+(count-1)*B_length:count*B_length);
end

%% pad with extra zeros to ensure that last gradient point is ZERO (sequence requirement)
g  = [afGradient; zeros(1,size(afGradient,2))];
rf = [rf; zeros(1,size(rf,2))];

%% 
iDimRF           = 2;                   % Amplitude/Phase of complex data
iDimGradient     = size(g,2);           % typically 3 (non-zero) gradient channels
iGradientSamples = length(g);
iSamples         = length(rf);

%% extract phase and amp of rf
abs_rfs = abs(rf);
angle_rfs = angle(rf);
index = angle_rfs<0;
angle_rfs(index) = angle_rfs(index)+2*pi;   % phase wrap

%% scaling and rounding 
g = g*10;  % [mT/m]
g         = round(g        *10^iFloatDigits)/10^iFloatDigits;
abs_rfs   = round(abs_rfs  *10^iFloatDigits)/10^iFloatDigits;
angle_rfs = round(angle_rfs*10^iFloatDigits)/10^iFloatDigits; 

%% find max rf and gradient
fMaxAbsRF = max(abs(rf(:)));   % maximum over all channels
fMaxAbsG1 = max(abs(g(:,1)));
fMaxAbsG2 = max(abs(g(:,2)));
fMaxAbsG3 = max(abs(g(:,3)));

% RFPulseID is contained in filename, typically starting from zero
strFileName = sprintf('%s%d.ini', PULSE_FILE_NAME, iRFPulseID);
fid = fopen( strFileName, 'w');

% header at start
fprintf(fid, '#pTXRFPulse - created by: <%s><%s>\n', FCTNAME, VERSION);

% here huge pseudocomment sections will be included 
%   for B0,B1 map, pulse design, SAR calculation ...

% global parameter section
fprintf(fid,'\n#----------\n');
fprintf(fid,'[pTXPulse]\n');
fprintf(fid,'\n');
fprintf(fid,'NUsedChannels    = %d\n', iNUsedChannels);
fprintf(fid,'DimRF            = %d\n', iDimRF);
fprintf(fid,'DimGradient      = %d\n', iDimGradient);
fprintf(fid,'MaxAbsRF         = %g\t\t # scaling for RF amplitude\n', fMaxAbsRF);
fprintf(fid,'\n');
fprintf(fid,'PulseName        = %s\t\t #standard RF pulse parameters\n', strPulseName);
fprintf(fid,'Family           = %s\n', strFamily);
fprintf(fid,'Comment          = %s\n', strComment);
fprintf(fid,'NominalFlipAngle = %g\n', fNominalFlipAngle);
fprintf(fid,'Samples          = %d\n', iSamples);
fprintf(fid,'AmplInt          = %g\n', fAmplInt);
fprintf(fid,'AbsInt           = %g\n', fAbsInt);
fprintf(fid,'PowerInt         = %g\n', fPowerInt);
fprintf(fid,'MinSlice         = %g\n', fMinSlice);
fprintf(fid,'MaxSlice         = %g\n', fMaxSlice);
fprintf(fid,'RefGrad          = %g\n', fRefGrad);





% gradient section
fprintf(fid,'\n#----------\n');
fprintf(fid,'[Gradient]\t # Gx\t Gy\t Gz\t [mT/m]\n');
fprintf(fid,'\n');
fprintf(fid,'GradientSamples   =  %g\n', iGradientSamples);
fprintf(fid,'MaxAbsGradient[0] =  %g\t %g\t %g\t\t # scaling for G amplitude\n', fMaxAbsG1, fMaxAbsG2, fMaxAbsG3);
fprintf(fid,'\n');
for iG = 1:iGradientSamples
    fprintf(fid,'G[%g]= %g\t %g\t %g\n', iG-1, g(iG,1), g(iG,2), g(iG,3));
end

% rf section
for iCh = 1:iNUsedChannels
    fprintf(fid,'\n#----------\n');
    fprintf(fid,'[pTXPulse_ch%g]\t # Amplitude\t Phase\n', iCh-1);
    fprintf(fid,'\n');
    for iRF = 1:iSamples
        fprintf(fid,'RF[%g]= %g\t %g\n', iRF-1, abs_rfs(iRF, iCh), angle_rfs(iRF, iCh));
    end
end
fprintf(fid,'\n');
fprintf(fid,'#EOF\n');
fclose(fid);

if( bVerbose )
    fprintf('%s(): <%s> successfully written. <%s>', FCTNAME, strFileName, VERSION);
end
bSuccess = 1;
return

%------------------------------------------------------------
%WARNING: need to make sure that gradient and rf go to zero at begin and
%end of the pulse, e.g. for spiral need to have rewinder
%
function bSuccess = create_header_parallel7T(g,b_est,oversample)

bSuccess = 0;          % global define

B_length = length(g)*oversample;
n_coil = length(b_est)/B_length;

for count = 1:n_coil
    rf(count,:) = b_est(1+(count-1)*B_length:count*B_length);
end


%DATA READING OUT TO FILES

max_x = max(g(:,1));
max_y = max(g(:,2));
max_z = max(g(:,3));
max_rf = max(abs(rf(:)));
RFpeakVolt = max_rf

min_x = min(g(:,1));
min_y = min(g(:,2));
min_z = min(g(:,3));
%min_rf = min(rf);

if abs(min_x)> abs(max_x)
    max_x = abs(min_x);
end
if max_x == 0
        max_x = 0.1;
end
if abs(min_y)> abs(max_y)
    max_y = abs(min_y);
end
if max_y == 0
        max_y = 0.1;
end
if abs(min_z)> abs(max_z)
    max_z = abs(min_z);    
end
if max_z == 0
        max_z = 0.1;
end
% if abs(min_rf)> abs(max_rf)
%     max_rf = abs(min_rf);    
% end
if max_rf == 0
    max_rf = 0.1;
end

angle_rfs = angle(rf);
index = angle_rfs<0;
angle_rfs(index) = angle_rfs(index)+2*pi;
% angle_rfs(8,:).'
% cumsum(g(:,3)*10)

for (ii = 1:8)
    
    % enter tx dir
    
    dirName = ['TX' num2str(ii)];
    eval(['cd ' dirName])
  
    % do grads
    fidx = fopen('gx.float', 'wb');
    fidy = fopen('gy.float', 'wb');
    fidz = fopen('gz.float', 'wb');

    fwrite(fidx, max_x*10, 'float32');
    fwrite(fidx, length(g(:,1)), 'float32');
    fwrite(fidx, g(:,1)*10, 'float32');

    fwrite(fidy, max_y*10, 'float32');
    fwrite(fidy, length(g(:,2)), 'float32');
    fwrite(fidy, g(:,2)*10, 'float32');

    fwrite(fidz, max_z*10, 'float32');
    fwrite(fidz, length(g(:,3)), 'float32');
    fwrite(fidz, g(:,3)*10, 'float32');

    fclose(fidx);
    fclose(fidy);
    fclose(fidz);

    % do rf amp phas
    
    fname_rfamp = sprintf('rf_amp.float');
    fidrf_amp = fopen(fname_rfamp, 'wb');
    fname_rfpha = sprintf('rf_pha.float');
    fidrf_pha = fopen(fname_rfpha, 'wb');
    fwrite(fidrf_amp, max_rf, 'float32');
    fwrite(fidrf_amp, length(rf), 'float32');
    fwrite(fidrf_amp, abs(rf(ii,1:end)), 'float32');
    fwrite(fidrf_pha, max_rf, 'float32');
    fwrite(fidrf_pha, length(rf), 'float32');
    fwrite(fidrf_pha, angle_rfs(ii,1:end), 'float32');
    fclose(fidrf_amp);
    fclose(fidrf_pha);
    
    cd ..
    
end
bSuccess = 1;
return

%--------------------------------------------------------------
function [dopt, narg] = handlerOptionsNarg(defaultOptions, inputOptions, narg, nargLimits)
% % 	x = x(dir, filenum, opt())
% % 	
% % 	handler for nargin and function options (optin)
% % 	
% % 	return: 
% % 	
% % 	Functions called: 
% % 	Tested for: 
% % 	
% %  	July 2004 -  Josef Pfeuffer
% % 	
% [dopt, narg] 	= handlerOptionsNarg(dopt, varargin, nargin, [0 nargVars]);

FCTNAME = 'handlerOptionsNarg';

      %%%%  option/narg handling   %%%%%
      % --- default options: define whole struct, comment
% defaultOptions

      % --- arg handling - error: narg > nargVars; demoMode: narg < nargVars
nargVars		= nargLimits(2);
error(nargchk(0,nargVars+1,narg)) 
if (narg == nargVars+1)
	if prod( size(inputOptions) ) > 1
		error(sprintf('%s: unexpected error<%s>', FCTNAME, size(inputOptions)))
	end
	optin	= inputOptions{1};
   	dopt 	= LocalSetopt(defaultOptions, optin);
   	narg 	= narg - 1;    % narg is NOT including options
else
	dopt 	= defaultOptions;
%	fprintf('%s: check nargLimits [%d %d]\n', FCTNAME, nargLimits(1), nargLimits(2))
end

return

%------------------------------------------------------------
function optout = LocalSetopt(optall, optin)
%
%%   updates default optall struct with input optin
%%   optout = setopt(optall, optin)
%%
%%   effects: warning on unkown option field
%%            ignores tagfield
%%
%%   JP Apr 2000

narg = nargin;
error(nargchk(2,2,narg));

tagfield = 'STRUCTNAME';
tagstr = 'optstruct';

optout = optall;
if (isempty(optin))
   return
end
descin = fieldnames(optin);
l_descin = length(descin);

for il=1:l_descin
   if ( strcmp(descin(il), tagfield) &  ...
	strcmp(getfield(optin,descin{il}), tagstr) )
      %% do nothing
      %% fprintf('field <%s> contains <%s>',descin{il}, ...
      %%             getfield(optin,descin{il}));
   elseif (  isfield(optall, descin(il)) )
      optout = setfield(optout, descin{il}, getfield(optin,descin{il}));
   else
      errmsg = sprintf('! option <%s> ignored !',descin{il});
      warning(errmsg);
   end
end % for

%--------------------------------------------------------------
function s = opt(varargin)
%
%%   builds struct from arguments
%%   optstruct = opt('OPTION1',var1,'OPTION2',var2, ...)
%%
%%   effects: fielddescriptors are converted to UPPERCASE
%%
%%   JP Apr 2000

tagfield = 'STRUCTNAME';
tagstr = 'optstruct';

len_arg = length(varargin);
s = [];
s = setfield(s,tagfield,tagstr);      % init s - return to call w/o argument
for ilen = 1:2:len_arg-1
   if (~ischar(varargin{ilen}))
      varargin
      error('opt(): descriptor is not a string !');
   end
   fielddescriptor = upper(num2str( varargin{ilen} ));
   s = setfield(s, fielddescriptor, varargin{ilen+1});
end

%------------------------------------------------------------
function bSuccess = demoMode(FCTNAME)
% %
fprintf('  %%\n  %%\n  %%   %s: entering DEMO mode\n  %%\n  %%\n',FCTNAME)

% demo data for testing
load PulseDesignData_090627
whos g rf

% fix wrong b_est size by zero padding
%if( rem(length(b_est),length(g_out)) ~= 0)     
%    b_est = [b_est; zeros( ceil(length(b_est)/length(g_out))*length(g_out)-length(b_est),1)];
%end

% create old TX format (TXi/*.float)    using 'opt('WRITEOLDFORMAT',1)'
% default: new INI format
bSuccess = create_pTXRFPulse( g, rf, opt('FACTOROVERSAMPLE',1,'NOMFLIPANGLE',90));

%------------------------------------------------------------
