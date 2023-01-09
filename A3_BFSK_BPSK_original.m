%** Combining 3 antennas, BFSK and BPSK**%
clear;
format long;
N = 8000;                                  %** Number of Transmitted symbols   **%
M = 16;
nbits = 4;                                 % Number of bits per symbol   %
met = zeros(1,16);
                                           % BPSK Bit   BFSK Bit  %
s1 = [1, 0];                               %  0           0         %
s2 = [0, 1];                               %  0           1         %
s3 = [-1, 0];                              %  1           0         %  
s4 = [0, -1];                              %  1           1         %
    
for l=1:11
%     np = 1.0*(l-1);                        %** Eb/N0 in dB
%     a = (0.5/nbits)*10^(-np/10);           %** Noise variance in linear                  **%
        np = 2.5*(l-1); 
     a = (0.5/4)*10^(-np/10); 
%     a=0;
%     h = crand(2,1);                        %** Fading channel (assumed Rayleigh fading)  **%
%     h = ones(2,1);                         %** Ideal channel  (AWGN only)                **%
%     h = crandn(4,1);
    es = 0;
    for k=1:1:N
      x=round(1*rand(1,4));                %** Input bits (3 bits per symbol)            **%
      h = crandn(4,1);
      if (x(3)==0 && x(4)==0)
          if (x(1) == 0 && x(2) == 0)   
              rd1 = h(1)*s1 + sqrt(a)*crandn(1,2);  %** received signal point for 0000  **%
          elseif (x(1) == 0 && x(2) == 1)
              rd1 = h(1)*s2 + sqrt(a)*crandn(1,2);  %** received signal point for 0100  **%
          elseif (x(1) == 1 && x(2) == 0)
              rd1 = h(1)*s3 + sqrt(a)*crandn(1,2);  %** received signal point for 1000  **%
          else
              rd1 = h(1)*s4 + sqrt(a)*crandn(1,2);  %** received signal point for 1100  **%
          end
          h(2)=h(1)*e^(j*pi/2)
      elseif (x(3)==1 && x(4)==0)
          if (x(1) == 0 && x(2) == 0)   
              rd1 = h(2)*s1 + sqrt(a)*crandn(1,2);  %** received signal point for 0001  **%
          elseif (x(1) == 0 && x(2) == 1)
              rd1 = h(2)*s2 + sqrt(a)*crandn(1,2);  %** received signal point for 0101  **%
          elseif (x(1) == 1 && x(2) == 0)
              rd1 = h(2)*s3 + sqrt(a)*crandn(1,2);  %** received signal point for 1001  **%
          else
              rd1 = h(2)*s4 + sqrt(a)*crandn(1,2);  %** received signal point for 1101  **%
          end
      elseif(x(3)==0 && x(4)==1)
           if (x(1) == 0 && x(2) == 0)   
              rd1 = h(3)*s1 + sqrt(a)*crandn(1,2);  %** received signal point for 0010  **%
          elseif (x(1) == 0 && x(2) == 1)
              rd1 = h(3)*s2 + sqrt(a)*crandn(1,2);  %** received signal point for 0110  **%
          elseif (x(1) == 1 && x(2) == 0)
              rd1 = h(3)*s3 + sqrt(a)*crandn(1,2);  %** received signal point for 1010  **%
          else
              rd1 = h(3)*s4 + sqrt(a)*crandn(1,2);  %** received signal point for 1110  **%
           end
      elseif(x(3)==1 && x(4)==1)
           if (x(1) == 0 && x(2) == 0)   
              rd1 = h(4)*s1 + sqrt(a)*crandn(1,2);  %** received signal point for 0011  **%
          elseif (x(1) == 0 && x(2) == 1)
              rd1 = h(4)*s2 + sqrt(a)*crandn(1,2);  %** received signal point for 0111  **%
          elseif (x(1) == 1 && x(2) == 0)
              rd1 = h(4)*s3 + sqrt(a)*crandn(1,2);  %** received signal point for 1011  **%
          else
              rd1 = h(4)*s4 + sqrt(a)*crandn(1,2);  %** received signal point for 1111  **%
           end
      end
    
      %**  Metric calculations    **%
      met(1) = sum(abs(rd1 - h(1)*s1).^2);      %  0 0 0 0%
      met(2) = sum(abs(rd1 - h(3)*s1).^2);      %  0 0 0 1%
      met(3) = sum(abs(rd1 - h(2)*s1).^2);      %  0 0 1 0%
      met(4) = sum(abs(rd1 - h(4)*s1).^2);      %  0 0 1 1%
      met(5) = sum(abs(rd1 - h(1)*s2).^2);      %  0 1 0 0%
      met(6) = sum(abs(rd1 - h(3)*s2).^2);      %  0 1 0 1%
      met(7) = sum(abs(rd1 - h(2)*s2).^2);      %  0 1 1 0%
      met(8) = sum(abs(rd1 - h(4)*s2).^2);      %  0 1 1 1%
      met(9) = sum(abs(rd1 - h(1)*s3).^2);      %  1 0 0 0%
      met(10) = sum(abs(rd1 - h(3)*s3).^2);      %  1 0 0 1%
      met(11) = sum(abs(rd1 - h(2)*s3).^2);      %  1 0 1 0%
      met(12) = sum(abs(rd1 - h(4)*s3).^2);      %  1 0 1 1%
      met(13) = sum(abs(rd1 - h(1)*s4).^2);      %  1 1 0 0%
      met(14) = sum(abs(rd1 - h(3)*s4).^2);      %  1 1 0 1%
      met(15) = sum(abs(rd1 - h(2)*s4).^2);      %  1 1 1 0%
      met(16) = sum(abs(rd1 - h(4)*s4).^2);      %  1 1 1 1%
      [med,nd] = min(met);

      out_sym = nd -1;                                               %** Output symbol decided by the decoder  **%
%       disp(out_sym);
      out_bit(1) = bitget(out_sym,4);                               %** First bit **%
      out_bit(2) = bitget(out_sym,3);                               %** Second bit  **%
      out_bit(3) = bitget(out_sym,2);
      out_bit(4) = bitget(out_sym,1);
      es = es + bitxor(x(1),out_bit(1)) + bitxor(x(2),out_bit(2)) + bitxor(x(3), out_bit(3) + bitxor(x(4), out_bit(4)));   %** Counting the number of wrong bits  **%
    end
    
%     ber = zeros(1,1000000);
%     snr = zeros(1,1000000);
    ber(l) = es/(4*N);                                                %** Calculating the average bit error probability  **%
    snr(l) = np;
    es = 0;
    srt = 10^(np/10);
     b_BFSK(l) = 0.5*erfc(sqrt(srt/2));
     b_BPSK(l) = 0.5*erfc(sqrt(srt));
    
    N = 200/ber(l);                                                   %** Number of transmitted symbols (updated)  **%
end

% semilogy(snr,ber,'bo-',snr,b_BFSK)
semilogy(snr,ber,'bo-', snr, b_BFSK, snr, b_BPSK)
grid
xlabel('SNR, in dB')
ylabel('BER')
legend('3-antenna SM','BFSK','BPSK')
axis([0 25 1e-6 1])