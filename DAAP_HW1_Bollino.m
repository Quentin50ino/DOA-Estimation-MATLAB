clear; close all; clc

% Parameters
M = 16;              % Number of microphones
L = 0.45;            % Array length in meters
c = 343;             % Speed of sound in m/s
d = L / (M - 1);     % Inter-element spacing
[y, Fs] = audioread("array_recordings.wav");

% STFT parameters
window_length = 256;
hop_size = 128;
nfft = 512;

theta = -90:1:90;

% Perform STFT
Y = stft_custom(y, window_length, hop_size, nfft);

% Initialize array to store pseudospectrum
P = zeros(length(theta), size(Y, 2));

% Process each time frame
for t = 1:size(Y, 2)
    % Process each frequency bin
    for f = 1:size(Y, 1)
        % Calculate steering vector
        freq = (f-1) * Fs / nfft;
        a = steering_vector(M, d, freq, theta, c);
        
        % Apply DAS beamformer
        P(:, t) = P(:, t) + abs(das_beamformer(squeeze(Y(f, t, :)), a)).^2;
    end
    
    % Normalize pseudospectrum (changed to logarithmic scale)
    P(:, t) = 10*log10(P(:, t) / max(P(:, t)));
end

% Visualization
visualize_pseudospectrum(P, theta, size(Y, 2), hop_size, Fs);
visualize_ula_and_doas(M, L, theta, P);
create_localization_video(M, L, theta, P, hop_size, Fs);

function Y = stft_custom(x, window_length, hop_size, nfft)
    % Custom Short-Time Fourier Transform implementation
    % x: input signal (can be multi-channel)
    % window_length: length of the analysis window
    % hop_size: hop size between consecutive frames
    % nfft: number of FFT points
    
    [num_samples, num_channels] = size(x);
    
    % Calculate number of time frames
    num_frames = floor((num_samples - window_length) / hop_size) + 1;
    
    % Initialize output
    Y = zeros(nfft/2+1, num_frames, num_channels);
    
    % Hann window
    win = hann(window_length);
    
    % Perform STFT
    for ch = 1:num_channels
        for frame = 1:num_frames
            start_idx = (frame-1)*hop_size + 1;
            end_idx = start_idx + window_length - 1;
            
            % Apply window and perform FFT
            frame_data = x(start_idx:end_idx, ch) .* win;
            spectrum = fft(frame_data, nfft);
            
            % Store only positive frequencies
            Y(:, frame, ch) = spectrum(1:nfft/2+1);
        end
    end
end

function a = steering_vector(M, d, freq, theta, c)
    % Calculate steering vector for ULA
    % M: number of microphones
    % d: inter-element spacing
    % freq: frequency
    % theta: array of angles (in degrees)
    % c: speed of sound
    
    k = 2 * pi * freq / c;  % Wave number
    m = (0:M-1)';  % Microphone indices
    
    a = exp(-1j * k * d * m * sind(theta));
end

function P = das_beamformer(Y, a)
    % Delay-and-Sum Beamformer
    % Y: STFT coefficients for one frequency bin and one time frame
    % a: steering vector
    
    P = abs(a' * Y).^2;
end

function visualize_pseudospectrum(P, theta, num_frames, hop_size, Fs)
    % Visualize averaged pseudospectrum
    figure;
    imagesc(theta, (0:num_frames-1) * hop_size / Fs, P');
    xlabel('DOA [deg]');
    ylabel('Time [s]');
    title('Frequency-averaged pseudospectrum');
    colorbar;
    set(gca, 'YDir', 'normal');
    
    % Adjust color scale
    colormap(jet);
    clim([-30 0]);  % Adjust this range as needed
    
    % Add gridlines
    grid on;
    
    % Adjust axis limits
    xlim([min(theta) max(theta)]);
    ylim([0 (num_frames-1) * hop_size / Fs]);
end

function visualize_ula_and_doas(M, L, theta, P)
    % Visualize ULA setup and estimated DOAs over time
    figure;
    
    % Plot ULA
    x = linspace(-L/2, L/2, M);
    y = zeros(1, M);
    plot(x, y, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    plot([-L/2, L/2], [0, 0], 'k-', 'LineWidth', 2);
    
    % Find DOA with maximum power for each time frame
    [~, max_idx] = max(P);
    doa = theta(max_idx);
    
    % Plot DOA arrows
    arrow_length = L/2;
    num_arrows = min(20, length(doa));  % Limit to 20 arrows or less
    step = floor(length(doa) / num_arrows);
    colors = jet(num_arrows);  % Create a colormap
    
    for i = 1:num_arrows
        t = (i-1)*step + 1;
        % Calculate arrow components
        dx = arrow_length * sind(doa(t));
        dy = arrow_length * cosd(doa(t));
        quiver(0, 0, dx, dy, 0, 'Color', colors(i,:), 'LineWidth', 2, 'MaxHeadSize', 0.5);
    end
    
    xlabel('x (m)');
    ylabel('y (m)');
    title('ULA Setup and Estimated Time-varying DOA');
    axis equal;
    grid on;
    xlim([-L/2, L/2]);  % Set x-axis limits to match the array length
    ylim([0, L/2]);  % Set y-axis limits symmetrically around the array
    
    % Add angle indicators
    angle_labels = [-90, -60, -30, 0, 30, 60, 90];
    for i = 1:length(angle_labels)
        angle = angle_labels(i);
        x_pos = (L/2) * sind(angle);
        y_pos = (L/2) * cosd(angle);
        plot([0, x_pos], [0, y_pos], 'k--', 'LineWidth', 1);
        text(x_pos, y_pos, sprintf('%d°', angle), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
    % Add colorbar to indicate time progression
    colormap(jet);
    c = colorbar;
    c.Label.String = 'Time progression';
    clim([1, num_arrows]);
end

function create_localization_video(M, L, theta, P, hop_size, Fs)
    % Create video showing DOA estimation over time
    fig = figure;
    
    % Find DOA with maximum power for each time frame
    [~, max_idx] = max(P);
    doa = theta(max_idx);
    
    % Create video writer object
    v = VideoWriter('localization_video.mp4', 'MPEG-4');
    v.FrameRate = 10;  % Adjust as needed
    open(v);
    
    for t = 1:size(P, 2)
        % Plot ULA
        x = linspace(-L/2, L/2, M);
        y = zeros(1, M);
        plot(x, y, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
        hold on;
        plot([-L/2, L/2], [0, 0], 'k-', 'LineWidth', 2);
        
        % Plot DOA arrow
        arrow_length = L/2;
        colors = jet(length(doa));  % Create a colormap
        dx = arrow_length * sind(doa(t));
        dy = arrow_length * cosd(doa(t));
        quiver(0, 0, dx, dy, 0, 'Color', colors(t,:), 'LineWidth', 2, 'MaxHeadSize', 0.5);
        
        xlabel('x (m)');
        ylabel('y (m)');
        title(sprintf('ULA Setup and Estimated DOA - Time: %.2f s', t * hop_size / Fs));
        axis equal;
        grid on;
        xlim([-L/2, L/2]);
        ylim([0, L/2]);
        
        frame = getframe(fig);
        writeVideo(v, frame);
        
        if(t ~= size(P, 2))
            clf;
        end
    end
    
    close(v);
end