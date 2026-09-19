%ΑΣΚΗΣΗ 4 ερώτημα 1
%pinakas A
A = [
    0 1 0 0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 1 0 1 0 1 0 0 0 0 0 0 0 0;
    0 1 0 0 0 1 0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 1 0 0 0;
    1 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 1 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 1 1 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 1 1 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 1 1 0 1 0 1;
    0 0 0 0 0 0 0 0 0 0 0 1 0 1 0
];


n = 15;
q = 0.15;
ni = sum(A, 2); 
G = zeros(n, n);
term2 = 0; 


for i = 1:n
    for j = 1:n

        if A(j, i) == 1
            term2 = (1 - q) / ni(j);
        else
            term2 = 0;
        end
        
        G(i, j) = (q / n) + term2;
        
    end
end

col_sums = sum(G, 1);
fprintf('\nΑθροίσματα στηλών πίνακα G\n');
disp(col_sums);

if all(abs(col_sums - 1) < 1e-10) && all(G(:) >= 0)
    fprintf('Άθροισμα στοιχείων κάθε στήλης = 1 και όλα τα στοιχεία του G >=0 \n');
else
    fprintf('\nΔεν είναι στοχαστικός\n');
end



%ΑΣΚΗΣΗ 4 --------------------------ερώτημα 2
p = ones(n, 1) / n; 

tol = 1e-12;
max_it = 1000;
error = 1;
iter = 0;
while error > tol && iter < max_it
    pPrevious = p;
    p = G * pPrevious;
    error = norm(p - pPrevious, inf); 
    iter = iter + 1;
end

for k = 1:n
    fprintf(' page %2d= %.10f\n', k, p(k));
end
fprintf('άθροισμα διανύσματος P= %.4f\n', sum(p));


% --- ASKHSH 4 erwthma 3


A2 = A;
A2(11, 1) = 1;
A2(13, 1) = 1;
A2(14, 1) = 1;
A2(15, 1) = 1;
A2(15, 12) = 0;

ni_2 = sum(A2, 2);
G2 = zeros(n, n);
term22 = 0;
for i = 1:n
    for j = 1:n
        
        if A2(j, i) == 1
            term22 = (1 - q) / ni_2(j);
        else 
            term22=0;
        end
        G2(i, j) = (q / n) + term22;
    end
end

p_2 = ones(n, 1) / n;
error = 1;
iter = 0;
while error > tol && iter < max_it
    pPrev = p_2;
    p_2 = G2 * pPrev;
    error = norm(p_2 - pPrev, inf);
    iter = iter + 1;
end

fprintf('PageRank σελίδας 1= %.10f\n', p(1));
fprintf('PageRank σελίδας 1= %.10f μετά από τις αλλαγές\n', p_2(1));

if p_2(1) > p(1)
    fprintf('-> βελτίωση σημαντικότητας\n');
end





%-------ΑΣΚΗΣΗ 4 ερώτημα 4-----------------
%pinakas A

A2 = A;
A2(11, 1) = 1; 
A2(13, 1) = 1; 
A2(14, 1) = 1; 
A2(15, 1) = 1; 
A2(15, 12) = 0;

n = 15;
q = 0.02;
ni = sum(A2, 2); 
G = zeros(n, n);
term2 = 0; 

for i = 1:n
    for j = 1:n

        if A2(j, i) == 1
            term2 = (1 - q) / ni(j);
        else
            term2 = 0;
        end
        
        G(i, j) = (q / n) + term2;
        
    end
end

col_sums = sum(G, 1);
fprintf('\nΑθροίσματα στηλών πίνακα G\n');
disp(col_sums);

if all(abs(col_sums - 1) < 1e-10) && all(G(:) >= 0)
    fprintf('Άθροισμα στοιχείων κάθε στήλης = 1 και όλα τα στοιχεία του G >=0 \n');
else
    fprintf('\nΔεν είναι στοχαστικός\n');
end

p = ones(n, 1) / n; 

tol = 1e-12;
max_it = 1000;
error = 1;
iter = 0;
while error > tol && iter < max_it
    pPrevious = p;
    p = G * pPrevious;
    error = norm(p - pPrevious, inf); 
    iter = iter + 1;
end

fprintf('\n q= 0.02 \n');
for k = 1:n
    fprintf(' page %2d= %.10f\n', k, p(k));
end
fprintf('άθροισμα διανύσματος P= %.4f\n', sum(p));


%--------------------------------------------------
fprintf('\n ---------------------- \n');

n = 15;
q = 0.6;
ni = sum(A2, 2); 
G = zeros(n, n);
term22 = 0; 

for i = 1:n
    for j = 1:n

        if A2(j, i) == 1
            term22 = (1 - q) / ni(j);
        else
            term22 = 0;
        end
        
        G(i, j) = (q / n) + term22;
        
    end
end

col_sums = sum(G, 1);
fprintf('\nΑθροίσματα στηλών πίνακα G\n');
disp(col_sums);

if all(abs(col_sums - 1) < 1e-10) && all(G(:) >= 0)
    fprintf('Άθροισμα στοιχείων κάθε στήλης = 1 και όλα τα στοιχεία του G >=0 \n');
else
    fprintf('\nΔεν είναι στοχαστικός\n');
end


p = ones(n, 1) / n; 

tol = 1e-12;
max_it = 1000;
error = 1;
iter = 0;
while error > tol && iter < max_it
    pPrevious = p;
    p = G * pPrevious;
    error = norm(p - pPrevious, inf); 
    iter = iter + 1;
end

fprintf('\n q= 0.6 \n');
for k = 1:n
    fprintf(' page %2d= %.10f\n', k, p(k));
end
fprintf('άθροισμα διανύσματος P= %.4f\n', sum(p));


%ΑΣΚΗΣΗ 4 ερώτημα 5-----------------------------
A2 = A;
A2(8, 11) = 3;
A2(12, 11) = 3;
disp(A2);

ni_2 = sum(A2, 2);
G2 = zeros(n, n);
q = 0.15;
term22 = 0;
for i = 1:n
    for j = 1:n
        
        if A2(j, i) > 0
            term22 = (A2(j,i) * (1 - q) )/ ni_2(j);
        else 
            term22=0;
        end
        G2(i, j) = (q / n) + term22;
    end
end

p1 = ones(n, 1) / n;
error = 1;
iter = 0;
while error > tol && iter < max_it
    pPrev = p1;
    p1 = G2 * pPrev;
    error = norm(p1 - pPrev, inf);
    iter = iter + 1;
end

fprintf('PageRank σελίδας 11= %.10f\n', p1(11));

% ΑΣΚΗΣΗ 4 ΕΡΩΤΗΜΑ 6
A3= A;
A3(10, :) = [];
A3(:, 10) = [];
n2 =14;
ni_3 = sum(A3, 2);
G3 = zeros(n2, n2);

for i = 1:n2
    for j = 1:n2
        if A3(j, i) == 1
            term3 = (1 - q) / ni_3(j);
        else
            term3 = 0; 
        end
        G3(i, j) = (q / n2) + term3;
    end
end

p2 = ones(n2, 1) / n2;
error = 1;
iter = 0;
while error > tol && iter < max_it
    pPrev = p2;
    p2 = G3 * pPrev;
    error = norm(p2 - pPrev, inf);
    iter = iter + 1;
end

for k = 1:n2
    fprintf(' page %2d= %.10f\n', k, p2(k));
end
fprintf('άθροισμα διανύσματος P= %.4f\n', sum(p2));


