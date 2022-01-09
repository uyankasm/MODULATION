    .................MUHAMMET ASIM UYANIK.................
......................1801022012.....................
clc;
clear all;

% Define Equations
Fs = 10^5;
Ts = 1/Fs;
t = 0:Ts:(5/50-Ts);
m = 20*cos(100*pi*t)+10*cos(200*pi*t); ...message signal
c = 100*cos(500*pi*t); ...carrier signal

% Message Signal m(t)
figure(1)
plot(t(1:2000),m(1:2000));
title('Message Signal m(t)')
xlabel('t(s)')
ylabel('|m(t)|')
grid on;

% Message Signal Fourier Transform m(f)
N = length(m);
f = linspace(-Fs/2-5,Fs/2-5,N);
MF = fftshift(fft(m)/N);
figure(2);
stem(f,abs(MF),'^');
title('Message Signal Fourier m(f)')
xlabel('f(hz)')
ylabel('|m(f)|')
grid on;
xlim([-200 200])

% Modulated Signal y(t)
yt = m.*c;
figure(3)
plot(t(1:2000),yt(1:2000));
title('Modulated Signal y(t)')
xlabel('t(sn)')
ylabel('|y(t)|')
grid on;

% y(t) Fourier Transform y(f)
Yf = fftshift(fft(yt)/N);
figure(4);
stem(f,abs(Yf),'^');
title('Modulated Signal Fourier y(f)')
xlabel('f(hz)')
ylabel('|y(f)|')
grid on;
xlim([-400 400])


% if 𝑐̂(𝑡) = 𝐶𝑜𝑠(500𝜋𝑡)...........................
c1 = cos(520*pi*t); ...𝑐̂(𝑡)

% e(t) before LPF
et = yt.*c1;
figure(5)
plot(t(1:2000),et(1:2000));
title('e(t) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(500𝜋𝑡) before LPF')
xlabel('t(sn)')
ylabel('|e(t)|')
grid on;

% e(t) Fourier Transform e(f)
Ef = fftshift(fft(et)/N);
figure(6);
stem(f,abs(Ef),'^');
title('e(f) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(500𝜋𝑡) before LPF')
xlabel('f(hz)')
ylabel('|e(f)|')
grid on;
xlim([-700 700]);

% LOW PASS FILTER...
A = zeros(1,49990/10);
A(1:10)=1;
B=zeros(1,50000/10);
B(end-9:end)=1;
C=[B 1 A];...FILTER

% z(f) low pass filter
Zf=Ef.*C;
figure(7);
stem(f,Zf,'^');
title('z(f) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(500𝜋𝑡) after LPF')
xlabel('f(hz)')
ylabel('|z(f)|')
grid on;
xlim([-200 200]);

% inverse fourier transform for z(f)
zt = N*ifft(ifftshift(Zf),N);
figure(8)
plot(t,zt);
title('z(t) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(500𝜋𝑡) after LPF')
xlabel('t(sn)')
ylabel('|z(t)|')
grid on;


% ÇALIŞMADI
% A=[f;abs(Ef)];
% [sat,sut]=size(A)
% 
% for i=1:sut+1
%     if A(1,i)<=-100 | A(1,i)>=100
%         A(2,i)=0;
%     end
% end
%         
% figure(7);
% stem(f,A(:,2));
% title('z(f) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(500𝜋𝑡) after LPF')
% xlabel('f(hz)')
% ylabel('|z(f)|')
% grid on;
% xlim([-1000 1000]);



% if 𝑐̂(𝑡) = 𝐶𝑜𝑠(520𝜋𝑡)..........................
c2 = cos(520*pi*t); ...new 𝑐̂(𝑡)

% e(t) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(520𝜋𝑡)  before LPF
et1 = yt.*c2;
figure(9)
plot(t(1:2000),et1(1:2000));
title('e(t) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(520𝜋𝑡) before LPF')
xlabel('t(sn)')
ylabel('|e(t)|')
grid on;

% e(t) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(520𝜋𝑡) Fourier Transform e(f)
Ef1 = fftshift(fft(et1)/N);
figure(10);
stem(f,abs(Ef1),'^');
title('e(f) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(520𝜋𝑡) before LPF')
xlabel('f(hz)')
ylabel('|e(f)|')
grid on;
xlim([-700 700])

% Low pass filter for e(f)
Zf1=Ef1.*C;
figure(11);
stem(f,Zf1,'^');
title('z(f) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(520𝜋𝑡) after LPF')
xlabel('f(hz)')
ylabel('|z(f)|')
grid on;
xlim([-200 200]);

% inverse fourier transform for z(f)
zt1 = N*ifft(ifftshift(Zf1),N);
figure(12)
plot(t,zt1);
title('z(t) @ 𝑐̂(𝑡) = 𝐶𝑜𝑠(520𝜋𝑡) after LPF')
xlabel('t(sn)')
ylabel('|z(t)|')
grid on;








