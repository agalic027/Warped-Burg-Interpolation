close all

[data, fs] = audioread('samples/audio_original.wav');
data = data(1 : round(length(data) / 3));
audio_spectrogram(data, fs, 'Original Signal')

damaged = [0.5; 0.55];
damaged_data = data;
damaged_data(round(fs * damaged(1)) : round(fs * damaged(2))) = 0;
audio_spectrogram(damaged_data, fs, 'Damaged Signal')

audiowrite('samples/audio_damaged.wav', damaged_data, fs);

model_order = 100;
lambda = 0.723;
alpha = 3;
interpolated_data = warped_burg_method ...
    (damaged_data, fs, damaged, model_order, lambda, alpha);
audio_spectrogram(interpolated_data, fs,['WBurg \lambda = ' num2str(lambda) ', p = ' num2str(model_order)])

audiowrite('samples/audio_interpolated.wav', interpolated_wburg, fs);
