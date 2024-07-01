clearvars;
pv = zeros(1,1);
X = load('network_example.txt');
pv(1,1) = RunPoisson(X);
% Print results
disp(['FCN p-value: ', num2str(pv(1,1))]);
