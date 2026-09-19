% -------------------ΑΣΚΗΣΗ 1 


% f και f'
f = @(x) 14*x.*exp(x-2) - 12*exp(x-2) - 7*x.^3 + 20*x.^2 - 26*x + 12;
df = @(x) (14*x + 2).*exp(x-2) - 21*x.^2 + 40*x - 26;
%γραφική παράσταση
x = linspace(0, 3, 1000);
figure;
plot(x, f(x), 'g');
hold on;
yline(0);
xlabel('x');
ylabel('f(x)');
legend('f(x)', 'y=0');

% Διχοτόμηση
fprintf('\n--- ΔΙΧΟΤΟΜΗΣΗ ---\n');
intervals = [0 0.5; 0.5 1; 1 1.5; 1.5 2; 1.2 2.5; 2.5 3; 0.8 1.5; 0 1.5; 0.5 1.5; 1.5 2.5];
orio = 0.5*1e-6;
for i = 1:10
    a = intervals(i, 1);
    b = intervals(i, 2);
    iters=0;
    if f(a) * f(b) < 0
        [root, iters] = bisection(f, a, b);
        fprintf('[%.1f, %.1f] : Ρίζα = %.6f,  Επαναλήψεις = %d\n', a, b, root, iters);
    elseif abs(f(a)) < orio 
        fprintf('[%.1f, %.1f] %.6f :ρίζα, επαναλήψεις = %d, πρόσημο θετικό ή 0\n', a, b, a, iters);
    elseif abs(f(b)) < orio 
        fprintf('[%.1f, %.1f] %.6f :ρίζα, επαναλήψεις = %d, πρόσημο θετικό ή 0 \n', a, b, b, iters);
    else
        fprintf('[%.1f, %.1f] η f διατηρεί πρόσημο  \n', a, b);
    end
end

% newton-raphson
fprintf('\n--- Newton-Raphson ---\n');
times = linspace(0, 3, 10);
for i = 1:10
    x0 = times(i);
    [root, iters] = newton_raphson(f, df, x0, 0.5*1e-6);
    fprintf('x0 = %.6f  →  Ρίζα = %.6f,  Επαναλήψεις = %d\n', x0, root, iters);
end

%τέμνουσας
fprintf('\n--- ΤΕΜΝΟΥΣΑ ---\n');
arxik = linspace (0, 3, 10);
for i = 1:10 
    x0=arxik(i);
    x1= x0 + 0.5;
    [root, iters] = secant(f, x0, x1);
    fprintf('x0 = %.6f, x1 = %.6f,  →  Ρίζα = %.6f,  Επαναλήψεις = %d\n', x0, x1, root, iters);
end

%--------------έλεγχος σύγκλισης----------------
[root1, iterations1, x_history1] = newton_raphson(f, df, 0.667, 0.5*1e-6);
[root2, iterations2, x_history2] = newton_raphson(f, df, 1.667, 0.5*1e-6);
foundRoots = [root1, root2];

for i = 1: length(foundRoots)
    r = foundRoots(i);
    fprintf('riza #%d: x = %.6f\n', i, r);
    fprintf('  f(%.6f)  = %.6e\n', r, f(r));
    fprintf('  f''(%.6f) = %.6e\n', r, df(r));
    
    if abs(df(r)) < 0.5*1e-6
        fprintf('(f'' κοντά στο 0):δεν γίνεται τετραγωνική σύγκλιση \n');
    else 
        fprintf('(f'' ≠ 0):γίνεται τετραγωνική σύγκλιση\n\n');
    end
end

error_ratios1 = [];
error_ratios2= [];

% riza 1
x0=0.5;
[root, iterations, x_history] = newton_raphson(f, df, x0, 0.5*1e-12);
fprintf('  riza : x = %.6f\n', root);
fprintf('  epanalipseis: %d\n', iterations);
errors= abs(x_history - root);
fprintf('  Sfalmata: %.6e\n', errors);

k=0;
    for j = 2 : length(errors)
        if errors(j-1)> 0.5*1e-10
            k=k+1;
            error_ratios1(k) = errors(j)/ (errors(j-1))^2 ;
        end
    end

%mesosOros
%tApoklish
mesosOros = mean(error_ratios1);
tApoklish = std(error_ratios1);
fprintf('mesos oros sfalmatwn gia riza 1: %.6e\n', mesosOros);
fprintf('tupikh apoklish sfalmatwn gia riza 1: %.6e\n', tApoklish);
if tApoklish < 0.5 && mesosOros >0 && mesosOros<100
    fprintf('STA8EROS LOGOS =  TETRAGWNIKH \n')
else
    fprintf('MH STA8EROS LOGOS =  OXI TETRAGWNIKH \n');
end
if tApoklish < 1.5 && mesosOros >0 && mesosOros<100
    fprintf('STA8EROS LOGOS =  TETRAGWNIKH  (me pio xalaro elegxo)\n')
else
    fprintf('MH STA8EROS LOGOS =  OXI TETRAGWNIKH \n');
end

%riza 2
x0=1.5;
[root, iterations, x_history] = newton_raphson(f, df, x0, 0.5*1e-12);
fprintf('  riza : x = %.6f\n', root);
fprintf('  epanalipseis: %d\n', iterations);
errors= abs(x_history - root);
fprintf('  sfalmata: %.6e\n', errors);

k=0;
    for j = 2 : length(errors)
        if errors(j-1)> 0.5*1e-10
            k=k+1;
            error_ratios2(k) = errors(j)/ (errors(j-1))^2 ;
        end
        
    end
mesosOros = mean(error_ratios2);
tApoklish = std(error_ratios2);
fprintf('mesos oros sfalmatwn gia riza 2: %.6e\n', mesosOros);
fprintf('tupikh apoklish sfalmatwn gia riza 2: %.6e\n', tApoklish);
if tApoklish < 0.5 && mesosOros >0 && mesosOros<100
    fprintf('STA8EROS LOGOS =  tetragnwikh \n');
else
    fprintf('MH STA8EROS LOGOS =  OXI TETRAGWNIKH \n');
end
if tApoklish < 1.5 && mesosOros >0 && mesosOros<100
    fprintf('STA8EROS LOGOS =  TETRAGWNIKH  (me pio xalaro elegxo)\n')
else
    fprintf('MH STA8EROS LOGOS =  OXI TETRAGWNIKH \n');
end






% διχοτόμησης
function [root, iterations] = bisection(f, a, b)
    maxIt = 1000;
    tol = 0.5*1e-6;
    if f(a) * f(b) > 0
        error('Δεν υπάρχει ρίζα στο [a,b]');
    end
    
    iterations = 0;
    while (b - a) / 2 > 0.5*1e-6 && iterations < maxIt
        m = (a+b)/2;  
        
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
    if iterations == maxIt
        root = NaN;
    else 
        root = (a+b) / 2;
    end  
end

%newton-raphson
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
            iterations = NaN;
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
    
    root = NaN;
    iterations = NaN;
end

% τέμνουσας
function [root, iterations] = secant(f, x0, x1)
    tol=0.5*1e-6;
    max_iter=100;
    
    iterations = 0;
    
    for i = 1:max_iter
        fx0 = f(x0);
        fx1 = f(x1);
        
       
        if abs(fx1) < 0.5*1e-12
            root = x1;
            return;
        end
        
        if abs(fx1 - fx0) < 0.5*1e-12
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

