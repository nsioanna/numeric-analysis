%-----------------άσκηση 2 ερώτημα 3


f = @(x) 54*x.^6 + 45*x.^5 - 102*x.^4 - 69*x.^3 + 35*x.^2 + 16*x - 4 ;
df = @(x) 324*x.^5 + 225*x.^4 - 408*x.^3 - 207*x.^2 + 70*x + 16 ;
ddf = @(x) 1620*x.^4 + 900*x.^3 - 1224*x.^2 - 414*x + 70 ;

x = linspace(-2, 2, 1000);
figure;
plot(x, f(x), 'g');
hold on;
yline(0);
xlabel('x');
ylabel('f(x)');
legend('f(x)', 'y=0');


fprintf(' --- NEWTON-RAPHSON ---');
% 100 τυχαίες αρχικές τιμές στο [-2, 2]
x0_list = -2 + 4*rand(100,1);

iters_classic = zeros(100,1);
iters_modified = zeros(100,1);

for i = 1:100
    x0 = x0_list(i);

    % Κλασική Newton
    [root1, it1] = newton_raphson(f, df, x0, 0.5*1e-7);
    if ~isnan(root1)
          iters_classic(i) = it1;
    else
          iters_classic(i) = NaN;
    end

    % Τροποποιημένη Newton
    [root2, it2] = newton_raphson2(f, df, ddf, x0, 0.5*1e-7);
    if ~isnan(root2)
          iters_modified(i) = it2;
    else
          iters_modified(i) = NaN;
    end
end

% Αποτελέσματα
mean_classic = mean(iters_classic,  'omitnan');
mean_modified = mean(iters_modified, 'omitnan');

std_classic = std(iters_classic , 'omitnan');
std_modified = std(iters_modified, 'omitnan');

fprintf('ΚΛΑΣΙΚΗ \n');
fprintf('ΜΕΣΟΣ ΟΡΟΣ ΕΠΑΝΑΛΗΨΕΩΝ : %.2f\n', mean_classic);
fprintf('ΤΥΠΙΚΗ ΑΠΟΚΛΙΣΗ : %.2f\n', std_classic);
fprintf(' ΤΡΟΠΟΠΟΙΗΜΕΝΗ \n');
fprintf('ΜΕΣΟΣ ΟΡΟΣ ΕΠΑΝΑΛΗΨΕΩΝ : %.2f\n', mean_modified);
fprintf('ΤΥΠΙΚΗ ΑΠΟΚΛΙΣΗ : %.2f\n', std_modified);


fprintf('\n --- ΔΙΧΟΤΟΜΗΣΗ --- \n');

brackets = [];
Ntarget = 100;

attempts = 0;
while size(brackets,1) < Ntarget && attempts < 20000
    attempts = attempts + 1;
    a = -2 + 4*rand();
    b = -2 + 4*rand();
    if a == b, continue; end
    if a > b, tmp = a; a = b; b = tmp; end
    fa = f(a); fb = f(b);
    if ~isfinite(fa) || ~isfinite(fb), continue; end
    if fa*fb < 0
        % optional: avoid duplicates (close intervals)
        too_close = false;
        for k=1:size(brackets,1)
            if abs(brackets(k,1)-a) < 1e-12 && abs(brackets(k,2)-b) < 0.5*1e-12
                too_close = true; break;
            end
        end
        if ~too_close
            brackets(end+1,:) = [a,b];
        end
    end
end

if size(brackets,1) < Ntarget
    warning('Δεν συγκεντρώθηκαν %d brackets, βρέθηκαν %d (attempts=%d)\n', Ntarget, size(brackets,1), attempts);
end

% shuffle order (προαιρετικό) για πιο "τυχαία" σειρά
idx = randperm(size(brackets,1));
brackets = brackets(idx,:);

numTests = min(100, size(brackets,1));
testIntervals = brackets(idx(1:numTests), :);

for i = 1:numTests
    a = testIntervals(i,1);
    b = testIntervals(i,2);

    [r, it] = bisection(f, a, b);  % η κλασική μέθοδος
    iters_dixotomhsh(i) = it;
end

iters_dixotomhsh2 = zeros(numTests,1);

for i = 1:numTests
    a = testIntervals(i,1);
    b = testIntervals(i,2);

    [r, it] = bisection2(f, a, b);
    iters_dixotomhsh2(i) = it;
end

mean_dixot = mean(iters_dixotomhsh, 'omitnan');
std_dixot  = std(iters_dixotomhsh, 'omitnan');

mean_dixot2 = mean(iters_dixotomhsh2, 'omitnan');
std_dixot2  = std(iters_dixotomhsh2, 'omitnan');

fprintf("ΚΛΑΣΙΚΗ ΔΙΧΟΤΟΜΗΣΗ\n");
fprintf("ΜΕΣΟΣ ΟΡΟΣ ΕΠΑΝΑΛΗΨΕΩΝ: %.2f\n", mean_dixot);
fprintf("ΤΥΠΙΚΗ ΑΠΟΚΛΙΣΗ: %.2f\n", std_dixot);

fprintf("ΤΡΟΠΟΠΟΙΗΜΕΝΗ ΔΙΧΟΤΟΜΗΣΗ\n");
fprintf("ΜΕΣΟΣ ΟΡΟΣ ΕΠΑΝΑΛΗΨΕΩΝ: %.2f\n", mean_dixot2);
fprintf("ΤΥΠΙΚΗ ΑΠΟΚΛΙΣΗ: %.2f\n", std_dixot2);



fprintf("\n--- ΤΕΜΝΟΥΣΑ ---\n");
% 100 τυχαίες αρχικές τιμές στο [-2, 2]
x0_list = -2 + 4*rand(100,1);
x1_list = -2 + 4*rand(100,1);
x2_list = -2 + 4*rand(100,1);

iters_secant  = zeros(100,1);
iters_secant2 = zeros(100,1);

for i = 1:100
    
    x0 = x0_list(i);
    x1 = x1_list(i);
    x2 = x2_list(i);

    % Κλασική τέμνουσα (χρησιμοποιεί μόνο x0, x1)
    [root1, it1] = secant(f, x0, x1);
    iters_secant(i) = it1;

    % Τροποποιημένη τέμνουσα (χρησιμοποιεί x0, x1, x2)
    [root2, it2] = secant2(f, x0, x1, x2);
    iters_secant2(i) = it2;

end

mean_secant  = mean(iters_secant, 'omitnan');
std_secant   = std(iters_secant, 'omitnan');

mean_secant2 = mean(iters_secant2, 'omitnan');
std_secant2  = std(iters_secant2, 'omitnan');



fprintf('--- ΚΛΑΣΙΚΗ ΤΕΜΝΟΥΣΑ ---\n');
fprintf('ΜΕΣΟΣ ΟΡΟΣ ΕΠΑΝΑΛΗΨΕΩΝ: %.2f\n', mean_secant);
fprintf('ΤΥΠΙΚΗ ΑΠΟΚΛΙΣΗ: %.2f\n', std_secant);

fprintf('--- ΤΡΟΠΟΠΟΙΗΜΕΝΗ ΤΕΜΝΟΥΣΑ ---\n');
fprintf("ΜΕΣΟΣ ΟΡΟΣ ΕΠΑΝΑΛΗΨΕΩΝ: %.2f\n", mean_secant2);
fprintf("ΤΥΠΙΚΗ ΑΠΟΚΛΙΣΗ: %.2f\n", std_secant2);





%διχοτόμηση klasikh
function [root, iterations] = bisection(f, a, b)
    max = 1000;
    tol = 0.5*1e-7;

    if f(a) * f(b) > 0
        root = NaN;
        iterations = NaN;
        return;	
    end
    
    iterations = 0;
    while (b - a) / 2 > tol && iterations < max
        m = (a+b)/2;  
        
        if f(m) == 0
            root = m; 
            return;
        end

        if abs(f(m)) < tol
            root = m;
            return;
        elseif f(a) * f(m) < 0 
            b = m; 
        else  
            a = m;  
        end
        
        iterations = iterations + 1;
    end
    
    root = (a+b) / 2;  
end

%διχοτόμηση tropopoihmenh
function [root, iterations] = bisection2(f, a, b)
    max = 1000;
    tol = 0.5*1e-7;

    if f(a) * f(b) > 0
        root = NaN;
        iterations = NaN;
        return;
    end

    iterations = 0;
    fa = f(a);
    fb = f(b);

    while (b-a)/2 > tol && iterations < max
        
        if abs(fa) < abs(fb)
            m = a + (b-a)/3 ;
        else    
            m = b - (b-a)/3 ; 
        end
        fm = f(m);
        if abs(fm) < tol
            root =m;
            return;
        end
        if fa * fm < 0  
            b = m;
            fb=fm;
        else                 
            a = m;
            fa = fm;
        end
    
        iterations = iterations+1;
    end
  
    root = (a+b) /2;
end


%newton-raphson klasikh
function [root, iterations, x_history] = newton_raphson(f, df, x0, tol) 
    max_iter=100;

    x = x0;
    iterations = 0;
    x_history = x0;
    
    for i = 1:max_iter
        fx = f(x);
        dfx = df(x);
        
        if abs(dfx) < 0.5*1e-12
            root = NaN;
            return;
        end
       
        x_new = x - fx/dfx;
        iterations = iterations + 1;
        x_history = [x_history; x_new];

        if abs(x_new - x) < tol
            root = x_new;
            return;
        end
        
        x = x_new;
    end
    
    root = NaN; % den brhke riza
    iterations = NaN;
end

%newton-raphson tropopoihmenh
function [root, iterations, x_history] = newton_raphson2(f, df, ddf, x0, tol) 
    max_iter=100;

    x = x0;
    iterations = 0;
    x_history = x0;
    
    for i = 1:max_iter
        fx = f(x);
        dfx = df(x);
        ddfx= ddf(x);
        
        if abs(fx) < 0.5*1e-12
            root = NaN;
            return;
        end   
        if abs(dfx) < 0.5*1e-12
            root = NaN;
            return;
        end
       
        a = (dfx/fx) - 0.5000000 * (ddfx/dfx) ;
        if abs(a) < 0.5*1e-12
            root = NaN;
            iterations = NaN;
            return;
        end
        x_new = x - 1/ a;  
        iterations = iterations + 1;
        x_history = [x_history; x_new];

        if abs(x_new - x) < tol
            root = x_new;
            return;
        end
        
        x = x_new;
    end
    
    root = NaN; % δεν βρήκε ρίζα
    iterations= NaN;
end

%Τέμνουσα κλασική
function [root, iterations] = secant(f, x0, x1)
    tol= 0.5*1e-7;
    max_iter=100;
    
    iterations = 0;
    
    for i = 1:max_iter
        fx0 = f(x0);
        fx1 = f(x1);
        
        % Έλεγχος αν η τιμή είναι μηδέν
        if abs(fx1) < tol
            root = x1;
            return;
        end
        
        if abs(fx1 - fx0) < tol
             root = NaN;
             iterations = NaN;
             return;
        end

        x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        iterations = iterations + 1;

        if abs(x_new - x1) < tol
            root = x_new;
            return;
        end
        
        x0 = x1;
        x1 = x_new;  
    end
    
    root = NaN;
    iterations = NaN;
end

%τέμνουσα τροποποιημένη
function [root, iterations] = secant2(f, x0, x1, x2)
    tol= 0.5*1e-7;
    max_iter=100;
    
    iterations = 0;
    
    for i = 1:max_iter
        fx0 = f(x0);
        fx1 = f(x1);
        fx2 = f(x2);
        
        if abs(fx0)< tol
            root = x0;
            return;
        end
        if abs(fx1)< tol
            root = x1;
            return;
        end
        if abs(fx2)< tol
            root = x2;
            return;
        end

        q = fx0 / fx1;
        r = fx2 / fx1;
        s = fx2 / fx0;
        paronomastis = (q-1)*(r-1)*(s-1) ; 

        if abs(paronomastis) < tol
            root = NaN;
            iterations = NaN;
            return;
        end

        arithmitis = r*(r-q)*(x2-x1) + (1-r)*s*(x2-x0); 
        x_new = x2 - arithmitis/paronomastis ;

        iterations = iterations + 1;

        if abs(x_new - x2) < tol
            root = x_new;
            return;
        end
        
        x0 = x1;
        x1 = x2;  
        x2 = x_new;  
    end
    
    root = NaN;
    iterations = NaN;
end














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
disp(max(max(abs(difference))));

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

% 1. από MATLAB
x_matlab = A \ b; 

% δικό μου
[P, L, U] = plu(A);
x_mine = solveLU(P, L, U, b);

% αποτελεσμ
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


% παράδειγμα για έλεγχο
n = 4;
B = rand(n);
A = B * B'; % Κόλπο για να φτιάξουμε σίγουρα θετικά ορισμένο πίνακα

% δικό μου
L_mine = cholesky(A);

% από MATLAB
L_matlab = chol(A, 'lower'); 

% πίνακες & σφάλματα
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