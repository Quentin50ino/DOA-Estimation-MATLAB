# DAAP Homework 1 – Source Localization with Delay-and-Sum Beamforming
## Purpose of the Project

This project aims to perform acoustic source localization using a Uniform Linear Array (ULA) with 16 microphones.
The implemented method is based on Delay-and-Sum (DAS) Beamforming, applied both in the time and frequency domains.

The work was developed as Homework 1 for the course Sound Analysis, Synthesis and Processing (DAAP).

## Implemented Features

ULA System: 16 microphones equally spaced over a 0.45 m array length.

Computation of inter-microphone spacing and use of sound speed c = 343m/s.

Sampling frequency: 8 kHz.

## Custom STFT:

Short-Time Fourier Transform implemented from scratch.

Parameters: Hann window (256 samples), hop size 128, FFT size 512.

## Beamforming (DAS):

Computation of steering vectors for angles from −90° to 90°.

Projection of microphone signals onto steering vectors.

Creation of the pseudospectrum in time and frequency.

## DOA Estimation (Direction of Arrival):

Identification of the arrival angle corresponding to the maximum power in the pseudospectrum.

## Visualizations:

Static plot of microphone array geometry and DOAs.

Video showing the dynamic evolution of source localization over time.

## Main Results

Frequency-averaged pseudospectrum illustrating power distribution over time and angle.

Dynamic DOA tracking, represented with arrows in static plots and animations.

Effective DAS beamforming demonstrated for real-time acoustic source localization.

Limitations and Considerations

Wideband sources: DAS is a narrowband technique → frequency averaging required for broadband signals.

Angular resolution: limited by the array length and microphone spacing.

## Conclusion

The project delivers a modular and robust framework for acoustic source localization using microphone arrays.
With the implementation of a custom STFT, DAS beamforming, and both static and dynamic visualizations, it provides a strong foundation for real-time sound source tracking applications.
