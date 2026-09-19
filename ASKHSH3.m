
%ΑΣΚΗΣΗ 3
%------------------άσκηση 3 ερώτημα 1

function [P, L, U] = plu(A)
    n = size(A, 1);
    L = eye(n);
    P = eye(n);  
    U = A;    
    
    for k = 1:n-1
        [maxVal, maxIdx] = max(abs(U(k:n, k)));
        pivot_row = maxIdx + k - 1;
        
        %πίνακας U
        temp = U(k, :);
        U(k, :) = U(pivot_row, :);
        U(pivot_row, :) = temp;
        
        %πίνακας P
        temp = P(k, :);
        P(k, :) = P(pivot_row, :);
        P(pivot_row, :) = temp;
        
        %πίνακας L
        if k > 1
            temp = L(k, 1:k-1);
            L(k, 1:k-1) = L(pivot_row, 1:k-1);
            L(pivot_row, 1:k-1) = temp;
        end
        
        % Gauss
        for i = k+1:n
            %συντελεστής
            factor = U(i, k) / U(k, k);
            L(i, k) = factor;
            %πίνακας U
            U(i, k:n) = U(i, k:n) - factor * U(k, k:n);
        end
    end
end



function x = solveLU(P, L, U, b)
    n = length(b);
    bNeo = P * b;
    y = zeros(n, 1);
    
    for i = 1:n
        y(i) = bNeo(i) - L(i, 1:i-1) * y(1:i-1);
    end
  
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end

%paradeigma
A = [2, -3, 1 ; 3, 1, -1 ; 1, -1, -1];
b = [1; 2; 1];

[P, L, U] = plu(A);
x = solveLU(P, L, U, b);

% Εκτύπωση αποτελεσμάτων
disp('Λύση x PRWTO:');
disp(x);

error = norm(A*x - b);
disp(['σφάλμα: ', num2str(error)]);


% παράδειγμα
A = [3 1 -4 1; 
     -5 2 1 -2; 
     -1 6 -3 -4; 
     -2 1 -4 2]; 

[P, L, U] = plu(A);

% έλεγχος
%PA - LU
difference = P*A - L*U;


disp('σφάλμα :');
disp(max(abs(difference)));

disp('P:'); disp(P);
disp('L:'); disp(L);
disp('U:'); disp(U);
b = [-3; 1; 2; 0];

x = solveLU(P, L, U, b);
disp('Λύση x:');
disp(x);
error = norm(A*x - b);
disp(['Σφάλμα: ', num2str(error)]);


A = [3 1 -4 1; 
     -5 2 1 -2; 
     -1 6 -3 -4; 
     -2 1 -4 2];

b = [-3; 1; 2; 0];

%από MATLAB
x_matlab = A \ b; 

%δικό μου
[P, L, U] = plu(A);
x_mine = solveLU(P, L, U, b);

%αποτελεσμ
disp('από MATLAB :');
disp(x_matlab);

disp('δικό μου:');
disp(x_mine);

disp('Διαφορά :');
disp(norm(x_matlab - x_mine));


%------------------άσκηση 3 ερώτημα 2

function L = cholesky(A)
    n = size(A, 1);
    L = zeros(n, n);
    
    %για κάθε στήλη
    for k = 1:n
        sum = 0;
        for j = 1:k-1
            sum = sum + L(k, j)^2;
        end
        
        val = A(k, k) - sum;
        if val <= 0
            error('μη θετικά ορισμένος πίνακας.');
        end
        L(k, k) = sqrt(val);
        
        for i = k+1:n
            sum2 = 0;
            for j = 1:k-1
                sum2 = sum2 + L(i, j) * L(k, j);
            end
            
            L(i, k) = (A(i, k) - sum2) / L(k, k);
        end
    end
end


%παράδειγμα για έλεγχο
n = 4;
B = rand(n);
A = B * B';

%δικό μου
L_mine = cholesky(A);

%από MATLAB
L_matlab = chol(A, 'lower'); 

%πίνακες & σφάλματα
disp('δικός μου:');
disp(L_mine);

disp(' A - L*L^T :');
disp(norm(A - L_mine * L_mine'));

disp('από MATLAB :');
disp(norm(L_mine - L_matlab));



%------------------άσκηση 3 ερώτημα 3


% gauss seidel
function x = func(n)
    x = zeros(n, 1);
    
    akriv = 0.5 * 1e-5;
    error = 1;
    
    iters = 0;
 
    while error > akriv
        iters = iters + 1;
        x_previous = x;
        
        %πρώτη γραμμή , β=3
        x(1) = (3 + 2*x(2)) / 5;
        
        %ενδιάμεσες γραμμές , β=1
        for i = 2:n-1
            x(i) = (1 + 2*x(i-1) + 2*x(i+1)) / 5;
        end
        
        %τελευταία γραμμή , β=3
        x(n) = (3 + 2*x(n-1)) / 5;
        
        %σφάλμα
        diff = abs(x - x_previous);
        error = max(diff);
    end
   
    fprintf('n=%d : %d επαναλήψεις.\n', n, iters);
end

%έλεγχος
fprintf(" n = 10 \n");
x10 = func(10);
disp(x10);

fprintf(" n = 5000 \n");
x5000 = func(5000);
fprintf(" πρώτες και τελευταίες 5 τιμές \n");
disp(x5000(1:5));
disp(x5000(4995:5000));