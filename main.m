bits_length = 1000;
Ts = 1;
Tb = 40;            %Ts = 1 & N = 40 --> Tb = N*Ts
N = 40;
No = 2;
Ws = (N*2*pi)/Tb;
bits = randi([0 1], 1, bits_length);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a random sequence of binary bits. -->PolarNRZ
figure;
[t, polar] = PolarNRZ(bits, Tb);
subplot(2, 1, 1);
plot(t, polar, 'LineWidth', 2.0);
axis([0 300 -3 3]);
grid on;
title('Polar NRZ Line Code');
subplot(2, 1, 2);
[PSDpolar, f] = periodogram(polar);   
plot(f, PSDpolar, 'LineWidth', 2.0);
axis([0 1 0 100]);
title("PSD of Polar NRZ Line Code");

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generate a random sequence of binary bits. -->UniPolarNRZ
figure;
[t, unipolar] = UniPolarNRZ(bits, Tb);
subplot(2, 1, 1);
plot(t, unipolar, 'LineWidth', 2.0);
axis([0 300 -3 3]);
grid on;
title('Uni-Polar NRZ Line Code');
subplot(2, 1, 2);
[PSDunipolar, f] = periodogram(unipolar);   
plot(f, PSDunipolar, 'LineWidth', 2.0);
axis([0 1 0 100]);
title("PSD of Uni-Polar NRZ Line Code");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a zero-mean, white Gaussian random sequence
noise = wgn(1, bits_length, 10*log10(No/2));
figure;
plot(noise, 'LineWidth', 2.0);
axis([0 300 -3 3]);
title('Noise');

w1 = 2*pi*(4+1)/Tb;
w2 = 2*pi*(4-1)/Tb;
wc = 8*pi/Tb;


t = Ts: Ts: Tb;

FSK_Error_Memeory = zeros(1, 9);
PSK_Error_Memeory = zeros(1, 9);

Theoretical_FSK_Error = zeros(1, 9);
Theoretical_PSK_Error = zeros(1, 9);

SNRe = -4:1:4;
T = Tb-t;
for SNR = 1:1:length(SNRe)
    A = sqrt((10^(SNRe(SNR)/10))*((2*No)/Tb));
    for iteration = 1:20
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate a sequence of FSK 
        FSK = [];
        for i = 1:1:length(bits)
            if bits(i) == 1
                Y = A*cos(w1*t);
            else
                Y = A*cos(w2*t);
            end
        FSK = [FSK Y];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate a sequence of PSK
        PSK = [];
        for i = 1:1:length(bits)
            if bits(i) == 1
                Y = A*cos(wc*t);
            else
                Y = -A*cos(wc*t);
            end
            PSK = [PSK Y];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting the PSK modulation result, the noise that will be
        % added and the addition result.
        noise = wgn(1, length(PSK), 10*log10(No/2));
        figure;
        subplot(3,1,1)
        plot(PSK);
        axis([0 300 -1 1]);
        title('Transmitted Signal - PSK');
        subplot(3,1,2)
        plot(noise);
        axis([0 300 -2 2]);
        title('Noise to be added');
        subplot(3,1,3)
        plot(PSK + noise);
        title('Received Signal - PSK + Noise');
        axis([0 300 -5 5]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting the FSK modulation result, the noise that will be
        % added and the addition result.
        noise = wgn(1, length(FSK), 10*log10(No/2));
        figure;
        subplot(3,1,1)
        plot(FSK);
        axis([0 300 -1 1]);
        title('Transmitted Signal - FSK');
        subplot(3,1,2)
        plot(noise);
        axis([0 300 -2 2]);
        title('Noise to be added');
        subplot(3,1,3)
        plot(FSK + noise);
        title('Received Signal - FSK + Noise');
        axis([0 300 -5 5]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Matched Filters
        h1_FSK = cos(w1*(Tb-t));
        h2_FSK = cos(w2*(Tb-t));
        h_PSK = cos(wc*(Tb-t));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Output of the matched-filter at the end of each bit period
        noise = wgn(1, length(FSK), 10*log10(No/2));
        MF_out_h1_FSK = conv((FSK + noise), h1_FSK);
        MF_out_h2_FSK = conv((FSK + noise), h2_FSK);
        
        noise = wgn(1, length(PSK), 10*log10(No/2));
        MF_out_h_PSK = conv((PSK + noise), h_PSK);

        MFS_out_h1_FSK = zeros(1, bits_length);
        MFS_out_h2_FSK = zeros(1, bits_length);
        MFS_out_h_PSK = zeros(1, bits_length);

%         for i = 1 : bits_length
%             MFS_out_h1_FSK(i*Tb) = MF_out_h1_FSK(i*Tb);
%             MFS_out_h2_FSK(i*Tb) = MF_out_h2_FSK(i*Tb);
%             MFS_out_h_PSK(i*Tb) = MF_out_h_PSK(i*Tb);
%         end

        for i = 1 : bits_length
            MFS_out_h1_FSK(i) = MF_out_h1_FSK(i*Tb);
            MFS_out_h2_FSK(i) = MF_out_h2_FSK(i*Tb);
            MFS_out_h_PSK(i) = MF_out_h_PSK(i*Tb);
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting the FSK matched filter output and the samples taken from
        % this output
        figure;
        subplot(2,2,1)
        plot(MF_out_h1_FSK);
        axis([0 300 -10 10]);
        title('MF FSK H1 Output Reciver');
        subplot(2,2,3)
        stem(MFS_out_h1_FSK);
        axis([0 300 -10 10]);
        title('MF FSK H1 Reciver Output Samples');
        subplot(2,2,2)
        plot(MF_out_h2_FSK);
        axis([0 300 -10 10]);
        title('MF FSK H2 Output Reciver');
        subplot(2,2,4)
        stem(MFS_out_h2_FSK);
        title('MF FSK H2 Reciver Output Samples');
        axis([0 300 -10 10]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting the PSK matched filter output and the samples taken from
        % this output
        figure;
        subplot(2,1,1)
        plot(MF_out_h_PSK);
        axis([0 300 -10 10]);
        title('MF PSK H Reciver Output');
        subplot(2,1,2)
        stem( MFS_out_h_PSK);
        axis([0 300 -10 10]);
        title('MF PSK H Reciver Output Samples');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FSK Decision Making
        FSK_Decision = zeros(1, bits_length);
        for i = 1:bits_length
            if (MFS_out_h1_FSK(i) > MFS_out_h2_FSK(i))
                FSK_Decision(i)=1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PSK Decision Making
        PSK_Decision = zeros(1, bits_length);
        for i = 1:bits_length
            if (MFS_out_h_PSK(i) > 0)
                PSK_Decision(i)=1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting the FSK Decision
        figure;
        subplot(2,1,1)
        stem(FSK_Decision);
        axis([0 300 -2 2]);
        title('FSK Decision');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting the PSK Decision
        subplot(2,1,2)
        stem(PSK_Decision);
        axis([0 300 -2 2]);
        title('PSK Decision');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate the bit-error rate for FSK
        total_number_of_erroneous_bits = 0;
        total_number_of_transmitted_bits = bits_length;
        for i = 1:total_number_of_transmitted_bits
            if bits(i) ~=  FSK_Decision(i)
                total_number_of_erroneous_bits = total_number_of_erroneous_bits + 1;
            end
        end
        FSK_Error_Memeory(SNR) = (total_number_of_erroneous_bits/total_number_of_transmitted_bits) + FSK_Error_Memeory(SNR);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate the bit-error rate for PSK
        total_number_of_erroneous_bits = 0;
        total_number_of_transmitted_bits = bits_length;
        for i = 1:total_number_of_transmitted_bits
            if bits(i) ~=  PSK_Decision(i)
                total_number_of_erroneous_bits = total_number_of_erroneous_bits + 1;
            end
        end
        PSK_Error_Memeory(SNR) = (total_number_of_erroneous_bits/total_number_of_transmitted_bits) + PSK_Error_Memeory(SNR);
        
%         X = sprintf('Using FSK and SNR = %d, the bit error is %d', SNRe(SNR), FSK_Error_Memeory(SNR));
%         disp(X);
%         X = sprintf('Using PSK and SNR = %d, the bit error is %d', SNRe(SNR), PSK_Error_Memeory(SNR));
%         disp(X);
        
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate the avg bit-error rate for FSK and PSK 
        FSK_Error_Memeory(SNR) = FSK_Error_Memeory(SNR)/20;
        PSK_Error_Memeory(SNR) = PSK_Error_Memeory(SNR)/20;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate the theoretical bit-error rate for FSK and PSK 
        Theoretical_FSK_Error(SNR) = 0.5 * erfc(sqrt(((A^2*Tb)/2)/No));
        Theoretical_PSK_Error(SNR) = 0.5 * erfc(sqrt(((A^2*Tb)/2)/(No*2)));
        
%         X = sprintf('Using FSK and SNR = %d, the avg bit error is %d', SNRe(SNR), FSK_Error_Memeory(SNR));
%         disp(X);
%         X = sprintf('Using PSK and SNR = %d, the avg bit error is %d', SNRe(SNR), PSK_Error_Memeory(SNR));
%         disp(X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the performance graphs
figure;
subplot(2,1,1);
semilogy(-4:1:4, FSK_Error_Memeory);
title('FSK Bit Error vs SNRe');
subplot(2,1,2);
semilogy(-4:1:4, PSK_Error_Memeory);
title('PSK Bit Error vs SNRe');

figure;
semilogy(-4:1:4, FSK_Error_Memeory);
hold on;
semilogy(-4:1:4, Theoretical_FSK_Error);
title('FSK Bit Error');
legend('Actual Bit Error','Theoretical Bit Error');
hold off;

figure;
semilogy(-4:1:4, PSK_Error_Memeory);
hold on;
semilogy(-4:1:4, Theoretical_PSK_Error);
title('PSK Bit Error');
legend('Actual Bit Error','Theoretical Bit Error');
hold off;


