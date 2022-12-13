%** Combined BFSK and BPSK Modulation                 **%

clear;
format long;
N = 8000;                                  %** Number of Transmitted symbols   **%
M = 8;
nbits = 3;                                 % Number of bits per symbol   %
met = zeros(1,8);
                                           % BPSK Bit   BFSK Bit  %
s1 = [1, 0];                               %  0           0         %
s2 = [0, 1];                               %  0           1         %
s3 = [-1, 0];                              %  1           0         %  
s4 = [0, -1];                              %  1           1         %

s5 = [1, 0];                               %  0           0        1 %
s6 = [0, 1];                               %  0           1        1 %
s7 = [-1, 0];                              %  1           0        1 %  
s8 = [0, -1];                              %  1           1        1 %
    
for l=1:11
    np = 1.0*(l-1);                        %** Eb/N0 in dB
    a = (0.5/nbits)*10^(-np/10);           %** Noise variance in linear                  **%
%     a=0;
%     h = crand(2,1);                        %** Fading channel (assumed Rayleigh fading)  **%
%     h = ones(2,1);                         %** Ideal channel  (AWGN only)                **%
    es = 0;
    for k=1:1:N
      x=round(1*rand(1,3));                %** Input bits (3 bits per symbol)            **%
      h = crandn(2,1);
      if (x(3)==0)
          if (x(1) == 0 && x(2) == 0)   
              rd1 = h(1)*s1 + sqrt(a)*crandn(1,2);  %** received signal point for 00              **%
          elseif (x(1) == 0 && x(2) == 1)
              rd1 = h(1)*s2 + sqrt(a)*crandn(1,2);  %** received signal point for 01              **%
          elseif (x(1) == 1 && x(2) == 0)
              rd1 = h(1)*s3 + sqrt(a)*crandn(1,2);  %** received signal point for 10              **%
          else
              rd1 = h(1)*s4 + sqrt(a)*crandn(1,2);  %** received signal point for 11              **%
          end
      else
          if (x(1) == 0 && x(2) == 0)   
              rd1 = h(2)*s1 + sqrt(a)*crandn(1,2);  %** received signal point for 00              **%
          elseif (x(1) == 0 && x(2) == 1)
              rd1 = h(2)*s2 + sqrt(a)*crandn(1,2);  %** received signal point for 01              **%
          elseif (x(1) == 1 && x(2) == 0)
              rd1 = h(2)*s3 + sqrt(a)*crandn(1,2);  %** received signal point for 10              **%
          else
              rd1 = h(2)*s4 + sqrt(a)*crandn(1,2);  %** received signal point for 11              **%
          end
      end
    
      %**  Metric calculations    **%
      met(1) = sum(abs(rd1 - h(1)*s1).^2);      %  0 0 0%
      met(2) = sum(abs(rd1 - h(2)*s1).^2);      %  0 0 1%
      met(3) = sum(abs(rd1 - h(1)*s2).^2);      %  0 1 0%
      met(4) = sum(abs(rd1 - h(2)*s2).^2);      %  0 1 1%
      met(5) = sum(abs(rd1 - h(1)*s3).^2);      %  1 0 0%
      met(6) = sum(abs(rd1 - h(2)*s3).^2);      %  1 0 1%
      met(7) = sum(abs(rd1 - h(1)*s4).^2);      %  1 1 0%
      met(8) = sum(abs(rd1 - h(2)*s4).^2);      %  1 1 1%
      [med,nd] = min(met);

      out_sym = nd -1;                                               %** Output symbol decided by the decoder  **%
%       disp(out_sym);
      out_bit(2) = bitget(out_sym,2);                                %** Second bit  **%
      out_bit(1) = bitget(out_sym,3);                                %** First bit **%
      out_bit(3) = bitget(out_sym,1);
      es = es + bitxor(x(1),out_bit(1)) + bitxor(x(2),out_bit(2)) + bitxor(x(3), out_bit(3));   %** Counting the number of wrong bits  **%
    end
    
%     ber = zeros(1,1000000);
%     snr = zeros(1,1000000);
    ber(l) = es/(3*N);                                                %** Calculating the average bit error probability  **%
    snr(l) = np;
    es = 0;
    srt = 10^(np/10);
     b_BFSK(l) = 0.5*erfc(sqrt(srt/2));
     b_BPSK(l) = 0.5*erfc(sqrt(srt));
    
    N = 200/ber(l);                                                   %** Number of transmitted symbols (updated)  **%
end

% semilogy(snr,ber,'bo-',snr,b_BFSK)
semilogy(snr,ber,'bo-', snr, b_BFSK)
grid
xlabel('SNR, in dB')
ylabel('BER')
axis([0 12 1e-10 1])
