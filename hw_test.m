% clc; clear all; close all;
hw02_worker = hw02();

x =sqrt(2);

[c, ~] = hw02_worker.p1(@p1_func, 1, 2, 1e-9, 'bisection');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p1_func, 1, 2, 1e-9, 'secant');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p1_func, 1, 2, 1e-9, 'newton', @p1_func_prime);
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p1_func, 1, 2, 1e-9, 'regula_falsi');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p1_func, 1, 2, 1e-9, 'steffensen');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p1_func, 1, 2, 1e-9, 'Illinois');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p1_func, 1, 2, 1e-9, 'Pegasus');
hw_assert(abs(c - x) < 1e-6);


x = 1.21;
[c, ~] = hw02_worker.p1(@p2_func, 1, 2, 1e-9, 'bisection');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p2_func, 1, 2, 1e-9, 'secant');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p2_func, 1, 2, 1e-9, 'newton', @p2_func_prime);
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p2_func, 1, 2, 1e-9, 'regula_falsi');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p2_func, 1, 2, 1e-9, 'steffensen');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p2_func, 1, 2, 1e-9, 'Illinois');
hw_assert(abs(c - x) < 1e-6);
[c, ~] = hw02_worker.p1(@p2_func, 1, 2, 1e-9, 'Pegasus');
hw_assert(abs(c - x) < 1e-6);


hw02_worker.p2();

function hw_assert(X)
    if X; fprintf('\t PASS\n'); else; fprintf('\t FAIL\n'); end
end

function ret = p1_func(x)
    ret = x^2 - 2;
end

function ret = p1_func_prime(x)
    ret = 2 * x;
end

function ret = p2_func(x)
    ret = sqrt(x) - 1.1;
end

function ret = p2_func_prime(x)
    ret = 0.5/sqrt(x);
end