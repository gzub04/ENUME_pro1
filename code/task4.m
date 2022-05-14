classdef task4
    properties
        A;          % initial matrix
        eigenvalues;% vector with eigenvalues
        n;          % size of matrix
        tolerance;  % upper bound for nulled elements
    end
    methods (Access = 'protected')

        % it's a function responsible for QR factorization taken from
        % the 68th page of the book Numerical Methods by Piotr Tatjewski
        % and slightly modified to fit into this class
        function [obj, Q, R]=qrmgs(obj)
            %QR (thin) factorization using modified Gram-Schmidt algorithm
            %for full rank real-valued and complex-valued matrices
            Q=zeros(obj.n,obj.n);
            R=zeros(obj.n, obj.n);
            d=zeros(1, obj.n);
            %factorization with orthogonal (not orthonormal) columns of Q:
            for i=1:obj.n
                Q(:,i)=obj.A(:,i);
                R(i,i)=1;
                d(i)=Q(:,i)'*Q(:,i);
                for j=i+1:obj.n
                    R(i,j)=(Q(:,i)'*obj.A(:,j))/d(i);
                    obj.A(:,j)=obj.A(:,j)-R(i,j)*Q(:,i);
                end
            end
            %column normalization (columns of Q orthonormal):
            for i=1:obj.n
                dd=norm(Q(:,i));
                Q(:,i)=Q(:,i)/dd;
                R(i,i:obj.n)=R(i,i:obj.n)*dd;
            end
        end

        % QR factorization without shift
        % page 77 from book Numerical Methods by Piotr Tatjewski
        function obj = eigval_QR_no_shift(obj)
            i = 1;
            loop_limit = 1000;
            while i <= loop_limit && max(max(obj.A - diag(diag(obj.A)))) > obj.tolerance    % loops until upper bound for nulled
                [obj, Q, R] = qrmgs(obj);   % qr factorization                              % elements isn't overstepped or lim reached
                obj.A = R * Q;              % store next iteration of A in obj.A; 
                i = i + 1;                  % in the end, obj.A will be a matrix with Eigenvalues on the diagonal
            end
            if i > loop_limit
                error('Loop limit reached. Program terminated.');
            else
                disp('Calculation successfull. Number of iterations:');
                disp(i);
            end
            obj.eigenvalues = diag(obj.A); % We store the eigenvalues in our eigenvector
        end

        % QR factorization with shift
        % page 77, Numerical Methods - Piotr Tatjewski
        % sligthly modified
        function obj = eigval_QR_w_shift(obj)
            INITIALsubmatrix = obj.A;   % initial matrix
            max_iter = 1000;        % max number of iteration for each eigenvalue
            total_iterations = 0;
            for k = obj.n:-1 : 2
                DK = INITIALsubmatrix; % DK - initial matrix to calculate a single eigenvalue
                i = 0;
                while i <= max_iter && max(abs(DK(k, 1:k-1))) > obj.tolerance   % loops until bottom row except diagonal value is = 0
                    DD = DK(k-1:k, k-1:k); % Bottom right 2x2 submatrix
                    tmp=[1,-(DD(1,1)+DD(2,2)),DD(2,2)*DD(1,1)-DD(2,1)*DD(1,2)]; % calculates roots 
                    eig_roots=roots(tmp);
                    if abs(eig_roots(1) - DD(2, 2)) < abs(eig_roots(2) - DD(2, 2))
                        shift = eig_roots(1);
                    else
                        shift = eig_roots(2);
                    end
                    DP = DK - eye(k) * shift;       % matrix shifted by pI
                    [Q1, R1] = external_qrmgs(DP);  % QR factorization on shifted matrix
                    DK = R1 * Q1 + eye(k) * shift;  % transformed matrix, shifted back
                    i = i + 1;
                end
                total_iterations = total_iterations + i; % update number of operations
                if i > max_iter
                    error('Loop limit exceeded program terminated');
                end
                obj.eigenvalues(k) = DK(k,k);       % save found eigenvalue
                if k > 2
                    INITIALsubmatrix = DK(1:k-1, 1:k-1);
                else
                    obj.eigenvalues(1) = DK(1, 1);  % Last eigenvalue is taken without looping again
                end
            end
            disp('Total number of iterations for all calculated matrices: ');
            disp(total_iterations);
        end
    end
    methods
        % constructor
        function obj = task4()
            obj.n = 5;
            obj.tolerance = 1e-6;
            obj.A = [10 18 -7 -37 2; 18 32 -2 3 14; -7 -2 13 -24 11; -37 3 -24 -37 21; 2 14 11 21 37]; % some symmetric matrix 5x5
            obj.eigenvalues = zeros(obj.n, 1);
        end

        % choose whether with or without shift
        function obj = find_eigenvalues(obj, shift)
            if shift == 'y'
                obj = eigval_QR_w_shift(obj);
            elseif shift == 'n'
                obj = eigval_QR_no_shift(obj);
            end
            disp('Eigenvalues:');
            disp(obj.eigenvalues);
        end
    end
end
