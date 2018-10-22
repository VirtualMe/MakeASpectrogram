% usage: [Timefreq [, f [, t]]] = specgram(x [, n [, Fs [, window [, overlap]]]])
% Original version by Paul Kienzle; modified by Sean Fulop March 2002

% Generate a spectrogram for the signal. This chops the signal into
% overlapping slices, windows each slice and applies a Fourier
% transform to determine the frequency components at that slice.
%
% x: vector of samples
% n: size of fourier transform window, or [] for default=256
% Fs: sample rate, or [] for default=2 Hz
% window: shape of the fourier transform window, or [] for default=hanning(n)
%    Note: window length can be specified instead, in which case
%    window=hanning(length)
% overlap: overlap with previous window, or [] for default=length(window)/2

% Example
%    x = chirp([0:0.001:2],0,2,500);  # freq. sweep from 0-500 over 2 sec.
%    Fs=1000;                  # sampled every 0.001 sec so rate is 1 kHz
%    step=ceil(20*Fs/1000);    # one spectral slice every 20 ms
%    window=ceil(100*Fs/1000); # 100 ms data window
%    specgram(x, 2^nextpow2(window), Fs, window, window-step);
%

function [Timefreq, f_r, t_r] = filespecgram(fname)
 if nargin < 1 
    [fname,pathname] = uigetfile('*.wav');
    fname = [pathname fname];
  end
[x, Fs] = audioread(fname);
x = x(:,1); % drop all but 1st channel
 n = 1024;

  % assign defaults
  window = hanning(n);
  overlap = length(window)/2;
  
  % if only the window length is given, generate hanning window
  if length(window) == 1, window = hanning(window); end

  % should be extended to accept a vector of frequencies at which to
  % evaluate the fourier transform (via filterbank or chirp
  % z-transform)
  if length(n)>1, 
    error('specgram doesn''t handle frequency vectors yet') 
end

  % compute window offsets
  win_size = length(window);
  if (win_size > n)
    n = win_size;
    warning('specgram fft size adjusted---must be at least as long as frame')
end
  step = win_size - overlap;

  % build matrix of windowed data slices
  offset = [ 1 : step : length(x)-win_size ];
  S = zeros (n, length(offset));
  for i=1:length(offset)
    S(1:win_size, i) = x(offset(i):offset(i)+win_size-1) .* window;
  end

  % compute fourier transform
  STFT = fft(S);
   
  % extract the positive frequency components
  if rem(n,2)==1
    ret_n = (n+1)/2;
  else
    ret_n = n/2;
end
  % inst = input('Do you want to see a spectrogram of this STFT (y/n)?', 's');
  inst = 'y';
  if strcmp(inst,'n') | strcmp(inst,'no')
      f = [0:ret_n-1]*Fs/n; t = offset/Fs;
    if nargout>0, Timefreq = STFT; end
    if nargout>1, f_r = f; end
    if nargout>2, t_r = t; end
 
      return
   
else 
     if nargout>0, Timefreq = STFT; end
 
  STFT = STFT(1:ret_n, :);

  f = [0:ret_n-1]*Fs/n;
  t = offset/Fs;
  
  if nargout>1, f_r = f; end
  if nargout>2, t_r = t; end

%  maxfreq = input('Please enter the maximum frequency (Hz)\n you wish to see on the display:');
% maxfreq = 10000;
  maxfreq = Fs/n; % start with df as max, go up as we find information
%  maglen = n*maxfreq/Fs;
  maglen = n/2;
  stftDim = size(STFT);
  if(stftDim(1)<maglen)
      maglen=length(STFT)
  end
  STFTmag = abs(STFT(2:maglen,:)); % magnitude of STFT
  % Aric - I added the nargin for maxp, and added the conditional here to
  % check whether it exists.  Original just did 
  % STFTmag = STFTmag/max(max(STFTmag)) for scaling.
  % max(max(STFTmag))
  % Aric - 7/18/18 - Change "maxp" to "refp" for it to have meaning when
  % using as a reference instead of normalization factor
  if nargin < 6 | isempty(refp)
	refp = 2e-5; % set normalization factor to Pa reference pressure
  end
  STFTmag = STFTmag/refp; % normalize so max magnitude will be 0 db
  STFTmag = max(STFTmag, 10^(-35/10)); % clip everything below -35 dB Still use -35 w/refp
    % attempt to adjust frequency axis on the fly - any frequency that goes
    % above 20 mPa or so is included...
    for fi=1:n/2-1
      if(max(STFTmag(fi,:))>100)
          maxfreq=fi*(Fs/n);
      end
    end

  
% imagesc(t, f(2:n*maxfreq/Fs), 20*log10(STFTmag)); axis xy; colormap(flipud(gray));
% display as an indexed grayscale image showing the log magnitude of the STFT, 
% i.e. a spectrogram; the colormap is flipped to invert the default setting,
% in which white is most intense and black least---in speech science we want
% the opposite of that.
[folder, filename, ext] = fileparts(fname);
hFig = figure('name',[filename ' Spectrogram'],'NumberTitle','off'); %added this so I could run consecutive plots Aric
subplot(5,1,1:4);
maxfi = floor(n*maxfreq/Fs);
pcolor(t, f(2:maxfi), 20*log10(STFTmag(1:maxfi-1,:)));
axis xy;
% Aric - I want the 2d colormap in color, not grayscale
colormap(jet);
% colormap(flipud(gray));
shading interp;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

sampletimes = 0:1/Fs:(size(x)-1)/Fs;
subplot(5,1,5);
plot(sampletimes,x,'DisplayName','Signal');
xlabel('Time (s)');
ylabel('Signal');

ButtonH=uicontrol('Parent',hFig, ...
    'Style','pushbutton',...
    'String','Play Audio Sc',...
    'Units','normalized',...
    'Position',[0.7 0.05 0.2 0.05],...
    'Visible','on',...
    'tag','playme',...
    'UserData',struct('data',x,'fs',Fs),...
    'Callback',@button_callback);

end 

return

end 

function button_callback(hObject,eventdata)
	h = findobj('tag','playme');
    soundsc(h.UserData.data,h.UserData.fs);
end