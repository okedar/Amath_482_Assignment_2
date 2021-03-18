% Clean workspace 
clear all; close all; clc

% In Order to run the program and extract calculations, this script must be
% run in sections such that only the data desired is loaded and only the
% computations desired are computed.

%%                                           % File Load In 
%%                                           
figure(1) % 14 seconds
[y, Fs] = audioread('GNR.m4a');
tr_gnr = length(y)/Fs; % record time in seconds 
plot((1:length(y))/Fs,y); 
xlabel('Time [sec]'); ylabel('Amplitude');
title('Sweet Child O Mine'); 
p8 = audioplayer(y,Fs); 
%playblocking(p8);

%%
figure(2) % 60 seconds
[y, Fs] = audioread('Floyd.m4a');
tr_gnr = length(y)/Fs; % record time in seconds 
plot((1:length(y))/Fs,y); 
xlabel('Time [sec]'); ylabel('Amplitude');
title('Comfortably Numb'); 
p8 = audioplayer(y,Fs); 
%playblocking(p8); 








%%                                            % Frequency and Time Domain song 1
%%
L = tr_gnr;   % time slot to transform
n = length(y); % number of Fourier modes 2^9
k = (1/L)*[0:n/2-1 -n/2:-1]; % Sweet child O Mine
k_shift = fftshift(k);
t2 = linspace(0,L,n+1); 
t = t2(1:n);
y_transpose = y.';
%%                                            % Ploting Frequencies song 1

%       Gabor and Spect
tau = 0:0.2:14; % centre of window
a = 200; % window size

for j = 1:length(tau)
    g_filter = exp(-a.*(t-tau(j)).^2); % Gaussian
    y_filtered = g_filter .* y_transpose;
    y_filt_trans = fft(y_filtered);
    Sgt_spec(:,j) = fftshift(abs(y_filt_trans));
end

figure(3)
pcolor(tau,k_shift,Sgt_spec)
shading interp
set(gca,'ylim',[210 760],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
yticks([277 311 370 415 554 698 739])
yticklabels({'C#', 'D#', 'F#', 'G#', 'C#', 'F', 'F#'})
title('Sweet Child O Mine'); 








%%                                            % Frequency and Time Domain song 2
%%
y_transpose = y.';
%y_transpose = y_transpose(1:1317960); % only for Comfortably Numb part 1
y_transpose = y_transpose(1317961:1317960*2); % only for Comfortably Numb part 2
tr_gnr = length(y_transpose)/Fs;
L = tr_gnr;   % time slot to transform
n = length(y_transpose); % number of Fourier modes 2^9
k = (1/L)*[0:n/2-1 -n/2:-1]; % Comfort Numb
k_shift = fftshift(k);
t2 = linspace(0,L,n+1); 
t = t2(1:n);
%%                                            % Ploting Frequencies song 2
tau = 0:0.2:30; % centre of window
a = 200; % window size

for j = 1:length(tau)
    g_filter = exp(-a.*(t-tau(j)).^2); % Gaussian
    y_filtered = g_filter .* y_transpose;
    y_filt_trans = fft(y_filtered);
    Sgt_spec(:,j) = fftshift(abs(y_filt_trans));
end

figure(3)
pcolor(tau,k_shift,Sgt_spec)
shading interp
set(gca,'ylim',[50 150],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('Comfortably Numb Part 2'); 
yticks([77 82 92 97.99 110 123])
yticklabels({'D#','E','G#','G','A','B'})
xticks([0 5 10 15 20 25 30]) % part 2
xticklabels({'30','35','40','45','50','55','60'}) % part 2








%%                                            % Frequency and Time Domain song 2 bass
%%
y_transpose = y.';
y_transpose = y_transpose(1:1317960); % only for Comfortably Numb
%y_transpose = y_transpose(1317961:1317960*2); % only for Comfortably Numb part 2
tr_gnr = length(y_transpose)/Fs;
L = tr_gnr;   % time slot to transform
n = length(y_transpose); % number of Fourier modes 2^9
k = (1/L)*[0:n/2-1 -n/2:-1]; % Sweet child O Mine
k = (1/L)*[0:n/2-1 -n/2:-1]; % Comfort Numb
k_shift = fftshift(k);
t2 = linspace(0,L,n+1); 
t = t2(1:n);
y_transform = fft(y_transpose);
shannon_filter = abs(k_shift) < 200;
y_trans_filt = y_transform .* fftshift(shannon_filter);
y_inverse = ifft(y_trans_filt);

    % Play Full Range then Guitar Only
y_guitar = y_inverse.';
%Full
p8 = audioplayer(y,Fs); 
playblocking(p8); 
%Guitar
p8 = audioplayer(y_guitar,Fs); 
playblocking(p8); 
%%                                            % Ploting Frequencies song 2 bass
tau = 0:0.2:30; % centre of window
a = 200; % window size

for j = 1:length(tau)
    g_filter = exp(-a.*(t-tau(j)).^2); % Gaussian
    y_filtered = g_filter .* y_inverse;
    y_filt_trans = fft(y_filtered);
    Sgt_spec(:,j) = fftshift(abs(y_filt_trans));
end

figure(3)
pcolor(tau,k_shift,Sgt_spec)
shading interp
set(gca,'ylim',[0 500],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('Comfortably Numb Part 1 Bass')








%%                                            % Frequency and Time Domain song 2 Guitar
%%
y_transpose = y.';
y_transpose = y_transpose(1:1317960); % only for Comfortably Numb
%y_transpose = y_transpose(1317961:1317960*2); % only for Comfortably Numb part 2
tr_gnr = length(y_transpose)/Fs;
L = tr_gnr;   % time slot to transform
n = length(y_transpose); % number of Fourier modes 2^9
k = (1/L)*[0:n/2-1 -n/2:-1]; % Sweet child O Mine
k = (1/L)*[0:n/2-1 -n/2:-1]; % Comfort Numb
k_shift = fftshift(k);
t2 = linspace(0,L,n+1); 
t = t2(1:n);
y_transform = fft(y_transpose);
shannon_filter = abs(k_shift) > 300;
shannon_filter1 = abs(k_shift) < 4000;
y_trans_filt = y_transform .* fftshift(shannon_filter).* fftshift(shannon_filter1);
y_inverse = ifft(y_trans_filt);

    % Play Full Range then Guitar Only
y_guitar = y_inverse.';
    %Full
%p8 = audioplayer(y,Fs); 
%playblocking(p8); 
    %Guitar
p8 = audioplayer(y_guitar,Fs); 
playblocking(p8); 
%%                                            % Ploting Frequencies song 2 Guitar
tau = 0:0.2:30; % centre of window
a = 200; % window size

for j = 1:length(tau)
    g_filter = exp(-a.*(t-tau(j)).^2); % Gaussian
    y_filtered = g_filter .* y_inverse;
    y_filt_trans = fft(y_filtered);
    Sgt_spec(:,j) = fftshift(abs(y_filt_trans));
end

figure(3)
pcolor(tau,k_shift,Sgt_spec)
shading interp
set(gca,'ylim',[0 5000],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('Comfortably Numb Part 1 Guitar')
% xticks([0 5 10 15 20 25 30]) % part 2
% xticklabels({'30','35','40','45','50','55','60'}) % part 2






