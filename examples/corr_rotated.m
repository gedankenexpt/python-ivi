function y = corr_rotated(downsampled, samples, theta, method)

rot = downsampled*exp(1j*(theta));
if (strcmpi(method,'cov'))
    % maximize covariance of I-data and Q-data
    y = 1-0.5*(mean(real(rot) .* real(samples)) + mean(imag(rot) .* imag(samples)));
    % or minimize covariance of I-Q and Q-I data
    %y = abs(mean(real(rot) .* imag(samples))) + abs(mean(imag(rot) .* real(samples)));
else
    evm = comm.EVM; %Read here https://se.mathworks.com/help/comm/ref/comm.evm-system-object.html;jsessionid=577b23876abe0c0d97d5e8687fff#bsnan5l-2_1 for more
    txdata=samples-mean(samples);
    rxdata = rot * rms(txdata)/rms(rot);
    rxdata = rxdata - mean(rxdata);
    y = step(evm,txdata,rxdata);    
end

end
