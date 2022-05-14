classdef task3    
properties
        A;                  % array A
        b;                  % array b
        n;                  % size of matrix
        x;                  % answers
        errors;             % errors in results
        iteration_errors;   % errors in equation after jaccobi/gauss-seidel method
        accuracy = 1e-10;   % accuracy required by the task
    end
    methods (Access = 'protected')
        % creates array for task 3 part 1
        function obj = task_3_create_array(obj)
            obj.A = [8 2 -3 1; 2 -25 5 -18; 1 3 15 -8; 1 1 -2 -10];
            obj.b = [7; 12; 24; 28];
        end

        % generates matrix A and b for task 2a
        function obj = task_a_create_arrays(obj)
            [row, col] = size(obj.A);
            % fill matrix A
            for i = 1 : row
                for j = 1 : col
                    if i == j
                        obj.A(i,j) = 9;
                    elseif (i == j - 1 || i == j + 1)
                        obj.A(i,j) = 1;
                    else
                        obj.A(i,j) = 0;
                    end
                end
                % fill matrix (vector) b
                obj.b(i) = 1.4 + 0.6 * i;
            end
        end

        % generates matrix A and b for task 2b
        function obj = task_b_create_arrays(obj)
            [row, col] = size(obj.A);
            % fill matrix A
            for i = 1 : row
                for j = 1 : col
                    obj.A(i, j) = 3/(4*(i + j - 1));
                end
                % fill matrix (vector) b
                if mod(i, 2) == 0       % if even
                    obj.b(i) = 0;
                else                    % if odd
                    obj.b(i) = 1/i;
                end
            end
        end

        % calculation error
        function obj = calculate_error(obj)
            for i = 1 : obj.n
                Ax = 0;        % sum to calculate Ax in a single row
                for j = 1 : obj.n
                    Ax = Ax + obj.A(i, j) * obj.x(j);
                end
                obj.errors(i) = Ax - obj.b(i);  % Ax - b saved in 'errors' vector
            end
        end

        % jacobi method
        function obj = jacobi(obj)
            error = inf;  % indicates current norm of error, set to inf so we enter the while loop at least once
            loop_limit = 1000;
            k = 1;          % number of iteration
            while (error > obj.accuracy && k <= loop_limit)   % loop until error is low enough or looped too many times
                x_prev = obj.x;         % we need to store previous x's as Jacobi's method uses only old iteration of x
                for i = 1 : obj.n       % i - row number
                    ax_sum = 0;
                    for j = 1 : obj.n   % summation sign in equation, j - column
                        if j ~=i        % we do not include diagonal values, as per definition
                            ax_sum = ax_sum + (obj.A(i,j) * x_prev(j));
                        end
                    end
                    obj.x(i) = (obj.b(i) - ax_sum) / obj.A(i, i);   % x_i = (b_i - sum(a_ij * x_j))/a_ii
                end
                obj = calculate_error(obj);     % calculate new errors
                error = norm(obj.errors);           % calculate the norm of the errors
                obj.iteration_errors(k) = error;    % and store them in vector
                k = k + 1;      % add to the iteration count
            end
        end

        % gauss-seidel method
        function obj = gauss_seidel(obj)
            error = inf;  % indicates current norm of error, set to inf so we enter the while loop at least once
            loop_limit = 1000;
            k = 1;          % number of iteration

            while (error > obj.accuracy && k <= loop_limit)   % loop until error is low enough or looped too many times
                for i = 1 : obj.n       % i - row number
                    ax_sum = 0;         % ax_sum will be the sum for the first and second summation
                    for j = 1 : obj.n   % summation signs in equation, j - column
                        if j ~=i        % we do not include diagonal values, as neither summation sign does
                            ax_sum = ax_sum + (obj.A(i,j) * obj.x(j));
                        end
                    end
                    obj.x(i) = (obj.b(i) - ax_sum) / obj.A(i, i);
                end
                obj = calculate_error(obj);     % calculate new errors
                error = norm(obj.errors);           % calculate the norm of the errors
                obj.iteration_errors(k) = error;    % and store them in vector
                k = k + 1;
            end
        end
    end
    methods
        % constructor
        function obj = task3(task_part)
            if task_part == 1       % depending on the part of the task, n differs
                obj.n = 4;
            elseif task_part == 2
                obj.n = 10;
            else
                disp('There is no such task part!');
            end

            if task_part == 1 || task_part == 2
                obj.A = zeros(obj.n, obj.n);    % the same as for task 2
                obj.b = zeros(obj.n, 1);
                obj.x = zeros(obj.n, 1);        % initial guess for iterations is 0 for all x
                obj.errors = zeros(obj.n, 1);
                obj.iteration_errors = [];
            end
        end

        % solves for task 3 part 1
        function obj = task3_p1_calculate(obj, method)
            obj = task_3_create_array(obj);     % creates array according to task's description
            if method == 'j'
                obj = jacobi(obj);
            elseif method == 's'
                obj = gauss_seidel(obj);
            else
                disp('Unknown solving method chosen in task3_calculate.');
            end
        end

        % solves for task 3 part 2
        function obj = task3_p2_calculate(obj, task_letter, method)
            if task_letter == 'a'   % first we choose task letter
                obj = task_a_create_arrays(obj);
            elseif task_letter == 'b'
                obj = task_b_create_arrays(obj);
            else
                disp('Wrong task letter in task_2_calculate!');
            end

            if (task_letter == 'a' || task_letter == 'b')
                if method == 'j'        % then we choose method that we want to use
                    obj = jacobi(obj);
                elseif method == 's'
                    obj = gauss_seidel(obj);
                else
                    disp('Uknown method chosen in task_2_calculate!');
                end
            end
        end        
    end
end
