%** Combined BFSK and BPSK Modulation                 **%

clear;
format long;
N = 8000;                                  %** Number of Transmitted symbols   **%
M = 4;
nbits = 2;                                 % Number of bits per symbol   %
                                           % BPSK Bit   BFSK Bit  %
s1 = [1, 0];                               %  0           0       %
s2 = [0, 1];                               %  0           1       %
s3 = [-1, 0];                              %  1           0       %  
s4 = [0, -1];                              %  1           1       %
    
for l=1:17
    np = 2.5*(l-1);                        %** Eb/N0 in dB
    a = (0.5/nbits)*10^(-np/10);           %** Noise variance in linear                  **%
    h = crand(2,1);                        %** Fading channel (assumed Rayleigh fading)  **%
    h = ones(2,1);                         %** Ideal channel  (AWGN only)                **%
    es = 0;
    for k=1:N
      x=round(1*rand(1,2));                %** Input bits (2 bits per symbol)            **%
	  if (x(1) == 0 && x(2) == 0)   
          rd1 = s1 + sqrt(a)*crandn(1,1);  %** received signal point for 00              **%
      elseif (x(1) == 0 && x(2) == 1)
          rd1 = s2 + sqrt(a)*crandn(1,1);  %** received signal point for 01              **%
      elseif (x(1) == 1 && x(2) == 0)
          rd1 = s3 + sqrt(a)*crandn(1,2);  %** received signal point for 10              **%
      else
          rd1 = s4 + sqrt(a)*crandn(1,2);  %** received signal point for 11              **%
      end
      %**  Metric calculations    **%
      met(1) = sum(abs(rd1 - h(1)*s1).^2);      %  0 0 %
      met(2) = sum(abs(rd1 - h(1)*s2).^2);      %  0 1 %
      met(3) = sum(abs(rd1 - h(1)*s3).^2);      %  1 0 %
      met(4) = sum(abs(rd1 - h(1)*s4).^2);      %  1 1 %
      [med,nd] = min(met);

      out_sym = nd -1;                                               %** Output symbol decided by the decoder  **%
      out_bit(2) = bitget(out_sym,1);                                %** First bit  **%
      out_bit(1) = bitget(out_sym,2);                                %** Second bit **%
      es = es + bitxor(x(1),out_bit(1)) + bitxor(x(2),out_bit(2));   %** Counting the number of wrong bits  **%
	end

    ber(l) = es/(2*N)                                                %** Calculating the average bit error probability  **%
    snr(l) = np;
    es = 0;
    srt = 10^(np/10);
    b_BFSK(l) = 0.5*erfc(sqrt(srt/2));
    b_BPSK(l) = 0.5*erfc(sqrt(srt));
    
    N = 200/ber(l)                                                   %** Number of transmitted symbols (updated)  **%
end

semilogy(snr,ber,'bo-',snr,b_BFSK)
grid
xlabel('SNR, in dB')
ylabel('BER')
legend('BPSK and BFSK','BFSK');
axis([0 40 1e-6 1])
