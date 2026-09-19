%ΑΣΚΗΣΗ 2
%------------------άσκηση 2 ερώτημα 1

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


% newton-raphson
fprintf('\n--- Newton-Raphson ---\n');
times = linspace(-2, 2, 10);
for i = 1:10
    x0 = times(i);
    [root, iters] = newton_raphson(f, df, ddf, x0, 0.5*1e-7);
    fprintf('%d ) x0 = %.7f  →  Ρίζα = %.7f,  Επαναλήψεις = %d\n',i, x0, root, iters);
end

%newton-raphson
function [root, iterations, x_history] = newton_raphson(f, df, ddf, x0, tol) 
    max_iter=100;

    x = x0;
    iterations = 0;
    x_history = x0;
    
    for i = 1:max_iter
        fx = f(x);
        dfx = df(x);
        ddfx= ddf(x);
        
        if abs(fx) < tol
            %fprintf('ρίζα = %.6f : , επαναλήψεις: %d \n', x, iterations);
            root = x;
            return;
        end   
        if abs(dfx) < tol
            root = NaN;
            iterations = NaN;
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
    
    root = NaN;
    iterations = NaN;
end

% DIXOTOMHSH
N = 5000;
xgrid = linspace(-2, 2, N);
brackets = [];
orio = 0.5*1e-7;

for i = 1:(N-1)
    a = xgrid(i);
    b = xgrid(i+1);
    fa = f(a);
    fb = f(b);

    if abs(fa) < orio
        brackets = [brackets; a a];
        continue;
    elseif abs(fb) < orio
        brackets = [brackets; b b];
        continue;
    end

    if fa * fb < 0
        brackets = [brackets; a b];
    end
end

fprintf("\nΒρέθηκαν %d πιθανά διαστήματα ριζών:\n", size(brackets,1));
for i = 1:size(brackets,1)
    if brackets(i,1) == brackets(i,2)
        fprintf('[%.7f] (Ακριβές σημείο)\n', brackets(i,1));
    else
        fprintf('[%.7f , %.7f] \n', brackets(i,1), brackets(i,2));
    end
end

intervals = brackets; 


fprintf("\nΔΙΧΟΤΟΜΗΣΗ\n");
for i = 1:size(intervals, 1)
    a = intervals(i, 1);
    b = intervals(i, 2);
    iters = 0;
    
    if a == b
        fprintf('Σημείο %.7f : Βρέθηκε ως ακριβής ρίζα από το grid search.\n', a);
        continue; 
    end

    if f(a) * f(b) <= 0
        [root, iters] = bisection(f, a, b); 
        fprintf('[%.7f, %.7f] : Ρίζα = %.7f,  Επαναλήψεις = %d\n', a, b, root, iters);
     elseif abs(f(a)) < orio 
        fprintf('[%.7f, %.7f] %.7f :ρίζα (στα άκρα), επαναλήψεις = %d\n', a, b, a, iters);
    elseif abs(f(b)) < orio 
        fprintf('[%.7f, %.7f] %.7f :ρίζα (στα άκρα), επαναλήψεις = %d\n', a, b, b, iters);
    else
        fprintf('[%.7f, %.7f] Προσοχή: Η f διατηρεί πρόσημο \n', a, b);
    end
end


function [root, iterations] = bisection(f, a, b)
    tol = 0.5*1e-7;

    if f(a) * f(b) > 0
        error('Δεν υπάρχει ρίζα στο [a,b]');
    end

    iterations = 0;

    while (b-a)/2 > tol
        fa = f(a);
        fb = f(b);
        
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
        else                 
            a = m;
        end
    
        iterations = iterations+1;
    end
  
    root = (a+b) /2;
end


%TEMNOUSA
%τέμνουσας ΚΛΑΣΙΚΗ -------------------
%fprintf('\n--- ΤΕΜΝΟΥΣΑ ---\n');
%arxik = linspace (-2, 2, 10);
%for i = 1:10 
%    x0=arxik(i);
%    x1= x0 + 0.5;
%    [root, iters] = secant(f, x0, x1);
%    fprintf('x0= %.7f, x1= %.7f,  →  Ρίζα = %.7f,  Επαναλήψεις = %d\n', x0, x1, root, iters);
%end

%τέμνουσας ΤΡΟΠΟΠΟΙΗΜΕΝΗ -------------------
fprintf('\n---ΤΕΜΝΟΥΣΑ ---\n');
arxik = linspace (-2, 2, 10);
for i = 1:10 
    x0=arxik(i);
    x1= x0 + 0.05;
    x2 = x1 + 0.05;
    [root, iters] = secant2(f, x0, x1, x2);
    fprintf('x0= %.7f, x1= %.7f, x2=%.7f  →  Ρίζα = %.7f,  Επαναλήψεις = %d\n', x0,x1,x2, root, iters);
end



% τέμνουσας
function [root, iterations] = secant2(f, x0, x1, x2)
    tol=0.5*1e-7;
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




% TEMNOUSA
function [root, iterations] = secant(f, x0, x1)
    tol=0.5*1e-7;
    max_iter=100;
    
    iterations = 0;
    
    for i = 1:max_iter
        fx0 = f(x0);
        fx1 = f(x1);
        
      
        if abs(fx1) < tol
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
