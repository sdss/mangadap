pro test_sn_definition

signal_1 = fltarr(100) + 10.   ;(constant signal of 10 counts, 100 pixels wide)
noise_1 = fltarr(100) + 3.   ;(constant error of 5 counts, 100 pixels wide)
signal_noise_1 = signal_1 / noise_1  ;(constant signal to noise vector (S/N per channel = 2), 100 pixels wide)

signal_2 = fltarr(1000) + 10.   ;(constant signal of 10 counts, 1000 pixels wide)
noise_2 = fltarr(1000) + 3.   ;(constant error of 5 counts, 1000 pixels wide)
signal_noise_2 = signal_2 / noise_2  ;(constant signal to noise vector (S/N per channel = 2), 1000 pixels wide)

; Ath formula: the mean signal to noise is the sum in quadrature of the
; signal to noise per channel, divided by the number of channels

mean_sn_1_A = total((signal_noise_1)^2) / n_elements(signal_noise_1)
mean_sn_2_A = total((signal_noise_2)^2) / n_elements(signal_noise_2)
print, mean_sn_1_A, mean_sn_2_A


;Bth formula: the mean signal to noise is the mean signal divided by
;the mean noise.

mean_signal_1 = total(signal_1) / n_elements(signal_noise_1)
mean_noise_1  =  sqrt( total(noise_1^2) / n_elements(noise_1) )
mean_sn_1_B = mean_signal_1/mean_noise_1

mean_signal_2 = total(signal_2) / n_elements(signal_noise_2)
mean_noise_2  =  sqrt( total(noise_2^2) / n_elements(noise_2) )
mean_sn_2_B = mean_signal_2/mean_noise_2

print,  mean_sn_1_B, mean_sn_2_B

stop
end
