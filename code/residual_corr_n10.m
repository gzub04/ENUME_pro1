function residual_corr_n10(task_letter, corr_count)
    matrix = task2(10);     % initialise matrix with n = 10
    if task_letter == 'a'
        matrix = task_a(matrix);
    elseif task_letter == 'b'
        matrix = task_b(matrix);
    else
        disp('Wrong letter!');
    end

    % display results and errors before corrections
    disp('Results and its errors before residual corrections:');
    for i = 1 : matrix.n
        fprintf("x%d = \t %f \t r%d = \t%f\n", i, matrix.x(i), i, matrix.errors(i));
    end

    %residual correction
    for i = 1 : corr_count
        matrix = residual_correction(matrix);
    end
    % print results after correction
    fprintf("Results and its errors after %d residual corrections:\n", corr_count);
    for i = 1 : matrix.n
        fprintf("x%d= \t %f \t r%d = %f\n", i, matrix.x(i), i, matrix.errors(i));
    end

    % condition number needed for explaining why matrix b is ill conditioned
    disp('Condition number:');
    disp(cond(matrix.A));
end
