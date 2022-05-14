classdef task2
    properties
        A;          % array A
        b;          % array b
        n;          % size of matrix
        x;          % answers
        errors;     % errors in results
    end
    methods (Access = 'protected')
        % generates matrix A and b for task a
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

        % generates matrix A and b for task b
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

        % gaussian elimination and partial pivoting
        function obj = gauss_and_partial_pivoting(obj)
            for i = 1 : obj.n                                   % main loop, going from the 1st to last row
                % pivoting is done first
                [max_num, max_row] = max(obj.A(i : obj.n, i));  % looks for the highest value in column i, below row i
                max_row = max_row + i - 1;  
                    % it is done to aquire number of row counting from the beginning, not from i, 
                    % as function max returns number of row starting from i
                if max_num == 0     % if column is all zeros, skip to the next iteration
                    continue;
                end
                if max_row ~= i     % if current row doesn't have highest value, swap:
                    obj.A([max_row i], :) = obj.A([i max_row], :);  % swaps row max_row with row i
                    obj.b([max_row i]) = obj.b([i max_row]);        % doesn't need a 2nd argument as it's a vector
                end
                % gaussian elimination
                for j = i + 1 : obj.n                   % loop for all rows, to find row multiplier and calculate rows using it
                    l = obj.A(j, i) / obj.A(i, i);      % calculate row multiplier
                    if l ~= 0                           % skip if multiplier is equal to 0
                        for k = i + 1 : obj.n           % loop to change each number in a single row using row multiplier
                            obj.A(j, k) = obj.A(j, k) - l * obj.A(i, k);
                        end
                        obj.b(j) = obj.b(j) - l * obj.b(j);
                    end
                end
            end
            % creating x vector with results
            obj.x(obj.n) = obj.b(obj.n) / obj.A(obj.n, obj.n);  % we start solving from the bottom; x(n) is trivial
            for i = obj.n - 1 : -1 : 1                          % i - rows; countdown from n - 1 to 1 with steps equal to -1
                sum = 0;
                for j = i + 1 : obj.n                           % j - column
                    sum = sum + obj.A(i, j) * obj.x(j);         % sum from the equation
                end
                obj.x(i) = (obj.b(i) - sum) / obj.A(i ,i);      % final result stored in x vector
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
    end
    methods
        % constructor, fills matrices with zeros
        function obj = task2(n)
            obj.A = zeros(n, n);
            obj.b = zeros(n, 1);
            obj.x = zeros(n, 1);
            obj.errors = zeros(n, 1);
            obj.n = n;
        end

        % creates matrix for task a
        function obj = task_a(obj)
            obj = task_a_create_arrays(obj);        % fills array according to task
            obj = gauss_and_partial_pivoting(obj);  % performs gaussian elimination with partial pivoting on the array
            obj = calculate_error(obj);             % calculates the error
        end

        % creates matrix for task b
        function obj = task_b(obj)
            obj = task_b_create_arrays(obj); % fills array according to task
            obj = gauss_and_partial_pivoting(obj);    % performs partial pivoting on the array
            obj = calculate_error(obj);     % calculates the error
        end
        
        % single residual correction
        function obj = residual_correction(obj)
            obj_copy = obj;
            obj_copy.b = obj_copy.errors;   % switch errors in b vector to fit into gauss_and_partial_pivoting(obj)
            obj_copy = gauss_and_partial_pivoting(obj_copy);    % solve the new obj
            obj.x = obj.x - obj_copy.x;     % new results according to x_2 = x_1 - delta(x)
            obj = calculate_error(obj);     % recalculate errors
        end
        

        % testing purposes, stores correct results (x) in 'errors' vector
        function obj = verify(obj)
            obj.errors = obj.A \ obj.b;
        end        
    end
end

