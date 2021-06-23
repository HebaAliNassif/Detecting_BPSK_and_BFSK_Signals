%Generating UniPolarNRZ line code  
function [t, x] = UniPolarNRZ(bits, Tb)
T = length(bits) * Tb;
N = 40;
Ws = (N*2*pi)/Tb;
Ts = (2*pi)/Ws;
t = 0:Ts:T;
x = zeros(1,length(t));
for i = 0:length(bits)-1
  if bits(i+1) == 1
    x(i*N+1:(i+1)*N) = 1;
  else
    x(i*N+1:(i+1)*N) = 0;
  end
end