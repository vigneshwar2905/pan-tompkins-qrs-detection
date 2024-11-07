# pan-tompkins-qrs-detection 
The Pan-Tompkins algorithm is widely used for QRS complex detection in ECG (electrocardiogram) signals. It involves a series of signal processing steps to enhance the QRS complex and suppress noise, allowing for precise identification of R-peaks. Below is a step-by-step guide on implementing this algorithm in MATLAB.

### Steps in the Pan-Tompkins Algorithm

1. **Bandpass Filtering:** Removes baseline wander and high-frequency noise.
2. **Differentiation:** Enhances the slope of the QRS complex.
3. **Squaring:** Emphasizes large differences and makes all data points positive.
4. **Moving Window Integration:** Smooths the squared signal to form peaks around QRS complexes.
5. **Thresholding and Peak Detection:** Detects R-peaks by applying a dynamic threshold.

### MATLAB Implementation

Here's how to implement each stage:

```matlab
% Load ECG data
fs = 360; % Sampling frequency, adjust based on your data
ecg_signal = load('ecg_data.mat'); % Load ECG data, modify with actual file name
ecg_signal = ecg_signal.ecg; % Assume 'ecg' is the variable name in the file

% 1. Bandpass Filter
% Using a Butterworth filter for bandpass (0.5 - 40 Hz)
low_cutoff = 0.5; % Low cutoff frequency
high_cutoff = 40; % High cutoff frequency
[b, a] = butter(1, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');
filtered_ecg = filtfilt(b, a, ecg_signal);

% 2. Differentiation
% Using a five-point derivative
diff_ecg = diff(filtered_ecg);
diff_ecg = [diff_ecg; 0]; % Append zero for equal length

% 3. Squaring
% Square each sample to amplify high-frequency components
squared_ecg = diff_ecg .^ 2;

% 4. Moving Window Integration
% Define window size, typically 150 ms, so window = 0.15 * fs
window_size = round(0.15 * fs);
integrated_ecg = movmean(squared_ecg, window_size);

% 5. Thresholding and Peak Detection
% Use a dynamic threshold for R-peak detection
threshold = 0.6 * max(integrated_ecg); % Initial threshold
qrs_peaks = find(integrated_ecg > threshold); % Identify peaks

% R-peak Refinement
r_peaks = [];
for i = 1:length(qrs_peaks)-1
    if qrs_peaks(i+1) - qrs_peaks(i) > round(0.2 * fs)
        [~, peak_loc] = max(filtered_ecg(qrs_peaks(i):qrs_peaks(i+1)));
        r_peaks = [r_peaks; qrs_peaks(i) + peak_loc - 1];
    end
end

% Plotting results
time = (0:length(ecg_signal)-1) / fs; % Time vector
figure;
plot(time, ecg_signal);
hold on;
plot(time(r_peaks), ecg_signal(r_peaks), 'ro');
title('ECG Signal with Detected R-Peaks');
xlabel('Time (s)');
ylabel('Amplitude');
legend('ECG Signal', 'Detected R-Peaks');
```

### Explanation of Code

1. **Bandpass Filter:** A Butterworth filter removes both low-frequency and high-frequency noise.
2. **Differentiation:** The derivative operation highlights the slopes of the QRS complex.
3. **Squaring:** The squaring operation amplifies high-frequency components and makes values positive.
4. **Moving Window Integration:** A sliding window smooths out the squared signal.
5. **Thresholding:** Peaks that exceed a certain threshold are identified as R-peaks. Adjust the threshold based on signal characteristics if necessary.

### Notes
- Adjust the `fs` (sampling frequency) according to your ECG data.
- Modify the `window_size` based on the sampling rate.
- Fine-tune the `threshold` variable to adapt to different ECG recordings.
