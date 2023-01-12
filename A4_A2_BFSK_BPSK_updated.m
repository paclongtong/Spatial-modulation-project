%** Combining 4 antennas, BFSK and BPSK**%
clear;
format long;
N = 100000;                                  %** Number of Transmitted symbols   **%
% N=1;
M = 16;

met = zeros(1,16);
                                           % BPSK Bit   BFSK Bit  %
s1 = [1, 0];                               %  0           0         %
s2 = [0, 1];                               %  0           1         %
s3 = [-1, 0];                              %  1           0         %  
s4 = [0, -1];                              %  1           1         %
    
nbits = 3;                                 % Number of bits per symbol   %
% ber2 = zeros(1,16);
for l=1:17
%     np = 1.0*(l-1);                        %** Eb/N0 in dB
%     a = (0.5/nbits)*10^(-np/10);           %** Noise variance in linear                  **%
    
    np = 2.5*(l-1);                               %** E_b/N_0 in dB = 10*log10(snr)  **%
    a = (0.5/3)*10^(-np/10);                      %** Variance of the thermal noise  N0/2 = 0.5(N0/Es) = 0.5(N0/Eb)/2 **%
%     a=0
    es = 0;
%     a=0;
%     h = crand(2,1);                        %** Fading channel (assumed Rayleigh fading)  **%
%     h = ones(2,1);                         %** Ideal channel  (AWGN only)                **%
%     es = 0;

    
    for k=1:1:N
      h = crandn(2,1);
      theta = angle(h);
      h(2) = h(2)*exp(1i*(theta(1)-theta(2)+pi/2));
      x=round(1*rand(1,3));                %** Input bits (3 bits per symbol)            **%
      
%       
%                   hold on 
%  quiver(0,0,real(h(1)),imag(h(1)));
%  quiver(0,0,real(h(2)),imag(h(2)));
%  text(real(h(1)),imag(h(1)),'h(1)');text(real(h(2)),imag(h(2)),'h(2)');
%  axis([-1,1,-1,1])
%  grid on

      

%             hold on 
%  quiver(0,0,real(h(1)),imag(h(1)));
%  quiver(0,0,real(h(2)),imag(h(2)));
%  text(real(h(1)),imag(h(1)),'h(1)');text(real(h(2)),imag(h(2)),'h(2)');
%  axis([-1,1,-1,1])
%  grid on
 
%       axis([-2,2,-2,2]);
%       plot(real(h(1)),imag(h(1)),'g*',real(h(2)),imag(h(2)),'bo');
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
      met2(1) = sum(abs(rd1 - h(1)*s1).^2);      %  0 0 0%
      met2(2) = sum(abs(rd1 - h(2)*s1).^2);      %  0 0 1%
      met2(3) = sum(abs(rd1 - h(1)*s2).^2);      %  0 1 0%
      met2(4) = sum(abs(rd1 - h(2)*s2).^2);      %  0 1 1%
      met2(5) = sum(abs(rd1 - h(1)*s3).^2);      %  1 0 0%
      met2(6) = sum(abs(rd1 - h(2)*s3).^2);      %  1 0 1%
      met2(7) = sum(abs(rd1 - h(1)*s4).^2);      %  1 1 0%
      met2(8) = sum(abs(rd1 - h(2)*s4).^2);      %  1 1 1%
      [med,nd] = min(met2);

      out_sym = nd -1;                                               %** Output symbol decided by the decoder  **%
%       disp(out_sym);
      out_bit(2) = bitget(out_sym,2);                                %** Second bit  **%
      out_bit(1) = bitget(out_sym,3);                                %** First bit **%
      out_bit(3) = bitget(out_sym,1);
      es = es + bitxor(x(1),out_bit(1)) + bitxor(x(2),out_bit(2)) + bitxor(x(3), out_bit(3));   %** Counting the number of wrong bits  **%
    end
    
%     ber = zeros(1,1000000);
%     snr = zeros(1,1000000);
    ber2(l) = es/(3*N);                                                %** Calculating the average bit error probability  **%
    snr(l) = np;
    es = 0;
%     srt = 10^(np/10);
%      b_BFSK(l) = 0.5*erfc(sqrt(srt/2));
%      b_BPSK(l) = 0.5*erfc(sqrt(srt));
    
%     if (ber2(l) ~= 0)
%         N = 200/ber2(l);                                                   %** Number of transmitted symbols (updated)  **%
%     else
%         if (l~=1)
%             N = 200/0.0001;
%         else
%             N = 8000;
%         end
%     end
    if(l >14)
        N=N*10;
    end
end

N = 100000;
nbits = 4;                                 % Number of bits per symbol   %
for l=1:17
%     np = 1.0*(l-1);                        %** Eb/N0 in dB
%     a = (0.5/nbits)*10^(-np/10);           %** Noise variance in linear                  **%
        np = 2.5*(l-1); 
     a = (0.5/nbits)*10^(-np/10); 
%      a=0;
%     a=0;
%     h = crand(2,1);                        %** Fading channel (assumed Rayleigh fading)  **%
%     h = ones(2,1);                         %** Ideal channel  (AWGN only)                **%
%     h = crandn(4,1);
    es = 0;
    
    for k=1:1:N
        h = crandn(4,1);
            theta = angle(h);
      h(2) = h(2)*exp(1i*(theta(1)-theta(2)+pi/2));
      h(3) = h(3)*exp(1i*(theta(1)-theta(3)+pi/4));
      h(4) = h(4)*exp(1i*(theta(1)-theta(4)-pi/4));

      x=round(1*rand(1,4));                %** Input bits (3 bits per symbol)            **%
      
      %%%%%%%%%%%%
%             hold on 
%  quiver(0,0,real(h(1)),imag(h(1)));
%  quiver(0,0,real(h(2)),imag(h(2)));
%  quiver(0,0,real(h(3)),imag(h(3)));
%  quiver(0,0,real(h(4)),imag(h(4)));
%  text(real(h(1)),imag(h(1)),'h(1)');text(real(h(2)),imag(h(2)),'h(2)');
%  text(real(h(3)),imag(h(3)),'h(3)');text(real(h(4)),imag(h(4)),'h(4)');
%  axis([-1,1,-1,1])
%  grid on

      
      
%       axis([-2,2,-2,2]);
%       plot(real(h(1)),imag(h(1)),'g*',real(h(2)),imag(h(2)),'bo',real(h(3)),imag(h(3)),'rd',real(h(4)),imag(h(4)),'bo');
      
%       hold on 
%  quiver(0,0,real(h(1)),imag(h(1)));
%  quiver(0,0,real(h(2)),imag(h(2)));
%  quiver(0,0,real(h(3)),imag(h(3)));
%  quiver(0,0,real(h(4)),imag(h(4)));
%  text(real(h(1)),imag(h(1)),'h(1)');text(real(h(2)),imag(h(2)),'h(2)');
%  text(real(h(3)),imag(h(3)),'h(3)');text(real(h(4)),imag(h(4)),'h(4)');
%  axis([-1,1,-1,1])
%  grid on
      
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
%      b_BFSK(l) = 0.5*erfc(sqrt(srt/2));
%      b_BPSK(l) = 0.5*erfc(sqrt(srt));
    
%     N = 200/ber(l);                                                   %** Number of transmitted symbols (updated)  **%
    
%     if (ber2(l) ~= 0)
%         N = 200/ber2(l);                                                   %** Number of transmitted symbols (updated)  **%
%     else
%         if (l~=1)
%             N = 200/ber2(l-1);
%         else
%             N = 8000;
%         end
%     end
    if(l >14)
        N=N*10;
    end
end

N=100000;
out_bit = zeros(2);
for l=1:17
%     np = 1.0*(l-1);                        %** Eb/N0 in dB
%     a = (0.5/nbits)*10^(-np/10);           %** Noise variance in linear                  **%
     np = 2.5*(l-1); 
     a = (0.5/2)*10^(-np/10); 
%      a=0;
%     h = crand(2,1);                        %** Fading channel (assumed Rayleigh fading)  **%
%     h = ones(2,1);                         %** Ideal channel  (AWGN only)                **%
    es = 0;
    
    
    for k=1:1:N
      h = crand(2,1);
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
      met1(1) = sum(abs(rd1 - h(1)*s1).^2);      %  0 0 %
      met1(2) = sum(abs(rd1 - h(1)*s2).^2);      %  0 1 %
      met1(3) = sum(abs(rd1 - h(1)*s3).^2);      %  1 0 %
      met1(4) = sum(abs(rd1 - h(1)*s4).^2);      %  1 1 %
      [med,nd] = min(met1);

      out_sym = nd -1;                                               %** Output symbol decided by the decoder  **%
      out_bit(2) = bitget(out_sym,1);                                %** First bit  **%
      out_bit(1) = bitget(out_sym,2);                                %** Second bit **%
      es = es + bitxor(x(1),out_bit(1)) + bitxor(x(2),out_bit(2));   %** Counting the number of wrong bits  **%
	end

    ber1(l) = es/(2*N);                                                %** Calculating the average bit error probability  **%
    snr(l) = np;
    es = 0;
    srt = 10^(np/10);
    b_BFSK(l) = 0.5*erfc(sqrt(srt/2))+0.5*(1-sqrt(srt/(2+srt)));
    b_BPSK(l) = 0.5*erfc(sqrt(srt))+0.5*(1-sqrt(srt/(1+srt)));
    
%     N = 200/ber1(l);                                                   %** Number of transmitted symbols (updated)  **%
%     if (ber2(l) ~= 0)
%         N = 200/ber2(l);                                                   %** Number of transmitted symbols (updated)  **%
%     else
%         if (l~=1)
%             N = 200/ber2(l-1);
%         else
%             N = 8000;
%         end
%     end
end


% semilogy(snr,ber,'bo-',snr,b_BFSK)
semilogy(snr,ber1,snr,ber2,snr,ber,'bo-', snr, b_BFSK, snr, b_BPSK)
legend('BPSK and BFSK','2-antenna SM','4-antenna SM', 'BFSK','BPSK');
grid
xlabel('SNR, in dB')
ylabel('BER')
title('Spatial modulation with optimized antenna model 1')
axis([0 40 1e-7 1])