clear;
format long;
N = 8000;                                        %** Number of transmitted bits   **%
M = 4;

s1 = 1;                                           %** Mapping of BPSK bit 1, unit energy  **%
s2 = -1;                                          %** Mapping of BPSK bit 0, unit energy  **%
    
for l=1:1:17
    np = 2.5*(l-1);                               %** E_b/N_0 in dB = 10*log10(snr)  **%
    a = (0.5/2)*10^(-np/10);                      %** Variance of the thermal noise  N0/2 = 0.5(N0/Es) = 0.5(N0/Eb)/2 **%
    es = 0;
    for k=1:1:N
      x=round(1*rand(1,2));                       %** Information Bits                  **%
      h = crandn(2,1);                            %** Fading Channel of antenna 1, 2    **%
	  if (x(1) == 0)                              %** x(1) is to select the antenna     **%
          if (x(2) == 0)                          %** x(2) is the BPSK bit              **%
            rd1 = h(1)*s1 + sqrt(a)*crandn(1,1);  %** rd1 is the received signal sample **%
          else
            rd1 = h(1)*s2 + sqrt(a)*crandn(1,1);
          end
      else
          if (x(2) == 0)   
            rd1 = h(2)*s1 + sqrt(a)*crandn(1,1);
          else
            rd1 = h(2)*s2 + sqrt(a)*crandn(1,1);
          end
      end
      met(1) = sum(abs(rd1 - h(1)*s1)^2);        %  0 0 %
      met(2) = sum(abs(rd1 - h(1)*s2)^2);        %  0 1 %
      met(3) = sum(abs(rd1 - h(2)*s1)^2);        %  1 0 %
      met(4) = sum(abs(rd1 - h(2)*s2)^2);        %  1 1 %
      [med,nd] = min(met);                       % Select the symbol that gives the minimum metric  **%

      out_sym = nd -1;
      out_bit(2) = bitget(out_sym,1);
      out_bit(1) = bitget(out_sym,2);
      es = es + bitxor(x(1),out_bit(1)) + bitxor(x(2),out_bit(2));
	end

    ber(l) = es/(2*N);                            %** Bit Error Probability Calculations  **%
    snr(l) = np;
    es = 0;
    srt = 10^(np/10);
    bes(l) = 0.5*erfc(sqrt(srt/2));
    N = 300/ber(l);                              %** Number of transmitted bits          **%
end

semilogy(snr,ber,'bo-')
grid
xlabel('$$E_b/N_0$$, in dB')
ylabel('Bit Error Probability')
axis([0 40 1e-6 1])
