using Plots

lowpass_1pole_fit_func(f, iR, f3dB) = iR ./ (1 .+ (f./f3dB).^2).^(1/2)

kB = 1.39e-23

johnsonCurrentNoise(T, R) = sqrt(4*kB*T*R)

temp = 3.54
resistance_ohm_sensor = 4
ir = johnsonCurrentNoise(temp, resistance_ohm_sensor)

cutoff_freq(R,L) = R/(2π*L)

resistance_ohm_filter = 1e3
inductance_filter = 1e-6
f3db = cutoff_freq(resistance_ohm_filter, inductance_filter)

# this was the measured value 
f3db = 20e3

xvals = range(1,10000,1000)
plot(xvals, lowpass_1pole_fit_func(xvals, ir, f3db), xscale=:log10)