function plot_errors(task_letter, n)
    % calculating number of plot points based on the equation n = 2^a * 10,
    % where a is number of points
    i = n/10;
    a = 1;
    while i ~= 1
        i = i/2;
        a = a + 1;
    end
    n_vector = zeros(a, 1);         % vector containing points for the graph
    error_norm_vector = zeros(a, 1);    % vector containing Euclidean norm
    for i = 1 : a
        n_vector(i) = 2^(i-1) * 10; % n = 10, 20, 40... can be rewritten as 2^a * 10
        matrix = task2(n_vector(i));% creates empty task2 object
        if task_letter == 'a'       % choose letter
            matrix = task_a(matrix);
        elseif task_letter == 'b'
            matrix = task_b(matrix);
        else
            disp('Wrong letter!');
        end
        error_norm_vector(i) = norm(matrix.errors); % Euclidean norm of errors is stored here
    end
    plot(n_vector, error_norm_vector);
end