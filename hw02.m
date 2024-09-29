% Author: Yevhenii Kovryzhenko / yzk0058@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw02

classdef hw02
    methods (Static)

        function [c, n] = p1(f, a, b, epsilon, name, f_prime)
            % p1: Implement numerical methods to find the root of a function
            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance
            % :param name: string, the name of the method
            % :param f_prime: function handle, the derivative of the function, only needed for Newton's method
            %
            % :return: 
            % c: real number, the root of the function
            % n: integer, the number of iterations
            if strcmp(name, 'bisection')
                n = 0;
                c = (a + b) / 2;
                while abs(f(c)) > epsilon
                    if f(a) * f(c) < 0
                        b = c;
                    else
                        a = c;
                    end
                    c = (a + b) / 2;
                    n = n + 1;
                end
            elseif strcmp(name, 'secant') % secant method
                % Write your code here
                n = 0;
                c = b;
                while abs(f(c)) > epsilon                    
                    c = c - (f(c)) * (c - a) / (f(c) - f(a));                    
                    a = b;
                    b = c;
                    n = n + 1;
                end
            elseif strcmp(name, 'newton') % Newton's method
                % Write your code here
                n = 0;
                c = a;
                while abs(f(c)) > epsilon                    
                    c = c - (f(c)) / f_prime(c);
                    n = n + 1;
                end
            elseif strcmp(name, 'regula_falsi') % false position method
                % Write your code here
                n = 0;
                c = (a + b) / 2;
                while abs(f(c)) > epsilon
                    if f(a) * f(c) < 0
                        b = c;
                    else
                        a = c;
                    end
                    c = a - f(a) * (b-a)/(f(b) - f(a));
                    n = n + 1;
                end
            elseif strcmp(name, 'steffensen') % Steffensen's method
                % Write your code here
                n = 0;
                c = a;
                g = @(x) (f(x + f(x)) - f(x)) / f(x);
                while abs(f(c)) > epsilon                    
                    c = c - (f(c)) / g(c);
                    n = n + 1;
                end
            elseif strcmp(name, 'Illinois')
                [c, n] = hw02.p3(f, a, b, epsilon);
            elseif strcmp(name, 'Pegasus')
                [c, n] = hw02.p4(f, a, b, epsilon);
            end
        end

        function p2()
            
        %     summarize the iteration number for each method name in the table
        
        %            name |  iter |      value |
        %       bisection |    32 |  +1.373471 |
        %          secant |     8 |  +1.373471 |
        %          newton |     5 |  +1.373471 |
        %    regula_falsi |    25 |  +1.373471 |
        %      steffensen |    13 |  +1.373471 |
        %        Illinois |     9 |  +1.373471 |
        %         Pegasus |     8 |  +1.373471 |

        methods = {'bisection', 'secant', 'newton', 'regula_falsi', 'steffensen', 'Illinois', 'Pegasus'};
        n_methods = length(methods);
        fprintf("%s %15s | %5s | %10s |\n", "%", "name", "iter", "value");

        f = @(x) 3*x^3 - 2*x^2 - 4;
        f_prime = @(x) 9*x^2 - 4*x;
        a = 1;
        b = 3;
        ftol = 1E-9;

        for i = 1:n_methods
            method = methods{i};
            [s, n] = hw02.p1(f, a, b, ftol, method, f_prime);
            fprintf("%s %15s | %5i | %+10f |\n", "%", method, n, s);
        end
            
        end

        function [x2, n] = p3(f, a, b, epsilon)
            % For 6630 only

            % Implement the Illinois method to find the root of a function

            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance

            % :return:
            % x2: real number, the root of the function
            % n: integer, the number of iterations

            % Write your code here

            x0 = a; x1 = b;
            f0 = f(x0); f1 = f(x1);
            n = 0;

            while true
                x2 = x0 - f0*(x1-x0)/(f1-f0);
                f2 = f(x2);
                n = n + 1;

                if abs(f2) < epsilon
                    return
                end

                while f1*f2 > 0
                    lambda = 1/2;
                    % x0 = x0; %why???
                    f0 = lambda*f0;

                    x1 = x2;
                    f1 = f2;
                    x2 = x0 - f0*(x1-x0)/(f1-f0);
                    f2 = f(x2);
                    n = n + 1;

                    if abs(f2) < epsilon
                        return
                    end
                end
                if f1*f2 < 0
                    x0 = x1;
                    f0 = f1;

                    x1 = x2;
                    f1 = f2;
                end

            end
        end

        function [x2, n] = p4(f, a, b, epsilon)
            % For 6630 only

            % Implement the Pegasus method to find the root of a function

            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance

            % :return:
            % x2: real number, the root of the function
            % n: integer, the number of iterations

            % Write your code here
            x0 = a; x1 = b;
            f0 = f(x0); f1 = f(x1);
            n = 0;

            while true
                x2 = x0 - f0*(x1-x0)/(f1-f0);
                f2 = f(x2);
                n = n + 1;

                if abs(f2) < epsilon
                    return
                end

                while f1*f2 > 0
                    lambda = f1 / (f1+f2);
                    % x0 = x0;
                    f0 = lambda*f0;

                    x1 = x2;
                    f1 = f2;
                    x2 = x0 - f0*(x1-x0)/(f1-f0);
                    f2 = f(x2);
                    n = n + 1;

                    if abs(f2) < epsilon
                        return
                    end
                end
                if f1*f2 < 0
                    x0 = x1;
                    f0 = f1;

                    x1 = x2;
                    f1 = f2;
                end

            end
        end
    end
end