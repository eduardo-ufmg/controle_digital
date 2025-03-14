clearvars; close all; clc

% Load setup parameters
setup_pH;

step = 10;
t0 = step;
tf_min = 50;
tf_seg = 60 * tf_min;
t = t0:step:tf_seg;

Ts = 40;
steps_in_sample = Ts/step;
T = t(1:steps_in_sample:end);

tstart = 600;

Q1 = 3 * ones(1, length(t));
Q3 = 2 * ones(1, length(t));

xc = repmat(x0, 1, length(t));
pHc = zeros(1, length(t));

xd = repmat(x0, 1, length(T));
pHd = zeros(1, length(T));

for k = 2:tstart/Ts+1
    kc = (k - 1) * steps_in_sample + 1;
    [xd(k), pHd(k), xc(kc:kc+steps_in_sample-1), pHc(kc:kc+steps_in_sample-1)] = simrk_pH(xc(kc-1), Q1(kc-1), Q3(kc-1), step, Ts, params, Kas);
end

k_test_start = k;

u2 = Q3;

du2 = zeros(size(u2));
du2(k_test_start:end) = 0.5 * ones(1, length(T) - k_test_start + 1);

u2 = u2 + du2;

for k = k_test_start+1:length(T)
    kc = (k - 1) * steps_in_sample + 1;
    [xd(k), pHd(k), xc(kc:kc+steps_in_sample-1), pHc(kc:kc+steps_in_sample-1)] = simrk_pH(xc(kc-1), Q1(kc-1), Q3(kc-1), step, Ts, params, Kas);
end

figure;
subplot(3,1,1);
plot(t, pHc, 'b', 'DisplayName', 'pHc');
hold on;
plot(T, pHd, 'r', 'LineWidth', 2, 'DisplayName', 'pHd');
xlabel('Time (s)');
ylabel('pH');
title('Simulated pH');
legend;
grid on;

subplot(3,1,3);
plot(t, u2, 'm', 'DisplayName', 'u2');
xlabel('Time (s)');
ylabel('Flow rate u2');
title('Flow rate u2');
legend;
grid on;
