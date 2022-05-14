function plot_task3(task_sign)
    if task_sign == 1
        % calculating for jacobi
        matrix = task3(1);      % creates empty matrix for task3
        matrix = task3_p1_calculate(matrix, 'j');   % fills and calculates matrix for jacobi method

        jacobi_errors = matrix.iteration_errors;    % saves errors from each iteration
        jacobi_no_of_iterations = length(matrix.iteration_errors);  % returns number of iterations

        disp('Results x for Jacobi:');
        disp(matrix.x);

        % calculating for gauss-seidel
        matrix = task3(1);
        matrix = task3_p1_calculate(matrix, 's');

        gauss_seidel_errors = matrix.iteration_errors;
        gauss_no_of_iterations = length(matrix.iteration_errors);
        
    elseif task_sign == 'a'
        matrix = task3(2);
        matrix = task3_p2_calculate(matrix, 'a', 'j');
        jacobi_errors = matrix.iteration_errors;
        jacobi_no_of_iterations = length(matrix.iteration_errors);

        disp('Results x for Jacobi:');
        disp(matrix.x);


        matrix = task3(2);
        matrix = task3_p2_calculate(matrix, 'a', 's');
        gauss_seidel_errors = matrix.iteration_errors;
        gauss_no_of_iterations = length(matrix.iteration_errors);
    elseif task_sign == 'b'
        matrix = task3(2);
        matrix = task3_p2_calculate(matrix, 'b', 'j');
        jacobi_errors = matrix.iteration_errors;
        jacobi_no_of_iterations = length(matrix.iteration_errors);

        disp('Results x for Jacobi:');
        disp(matrix.x);


        matrix = task3(2);
        matrix = task3_p2_calculate(matrix, 'b', 's');
        gauss_seidel_errors = matrix.iteration_errors;
        gauss_no_of_iterations = length(matrix.iteration_errors);
    else
        disp('Wrong task sign!');
    end
    if task_sign == 1 || task_sign == 'a' || task_sign == 'b'
        
        fprintf("Number of operations required for Jacobi method: %d\n\n\n", jacobi_no_of_iterations);
        disp('Results x for Gauss-Seidel:');
        disp(matrix.x);
        fprintf("Number of operations required for Gauss-Seidel method: %d\n", gauss_no_of_iterations);

        plot(1:jacobi_no_of_iterations, jacobi_errors, 1:gauss_no_of_iterations, gauss_seidel_errors);
        legend('Jacobi method', 'Gauss-Seidel method');
        xlabel('Number of iterations');
        ylabel('Norm');
    end
end
