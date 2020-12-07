clear all; close all; clc

k = linspace(0,10,1001);

FF = zeros(1,length(k));
GG = zeros(1,length(k));
for i = 1:length(k)
    [FF(i), GG(i)] = NewTheodorsenFunction(k(i));
end

figure(1)
plot(k,FF);
title('Real Part of C(k) vs. Reduced Frequency');
xlabel('Reduced Frequency (k)');
ylabel('Real Part of Theodorsen Function');

figure(2)
plot(k,GG)
title('Imaginary Part of C(k) vs. Reduced Frequency');
xlabel('Reduced Frequency (k)');
ylabel('Imaginary Part of Theodorsen Function');

figure(3)
plot(FF,GG)
title('Imaginary vs. Real Part of Theodorsen Function');
xlabel('Real Part of C(k)');
ylabel('Imaginary Part of C(k)');