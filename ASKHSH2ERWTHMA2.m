%------------------άσκηση 2 ερώτημα 2
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


fprintf('\n --- ΤΡΟΠΟΠΟΙΗΜΕΝΗ ΔΙΧΟΤΟΜΗΣΗ --- \n');

brackets = [];
Ntarget = 20;


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
            if abs(brackets(k,1)-a) < 1e-12 && abs(brackets(k,2)-b) < 1e-12
                too_close = true; break;
            end
        end
        if ~too_close
            brackets(end+1,:) = [a,b];
        end

        len = abs(a-b);
        [root, iters] = bisection2(f, a, b);
        fprintf('[%.7f, %.7f] : Ρίζα = %.7f,  Επαναλήψεις = %d, μήκος = %.7f\n', a, b, root, iters, len);
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

iters_dixotomhsh2 = zeros(numTests,1);

for i = 1:numTests
    a = testIntervals(i,1);
    b = testIntervals(i,2);

    [r, it] = bisection2(f, a, b);
    iters_dixotomhsh2(i) = it;
end

mean_dixot2 = mean(iters_dixotomhsh2 , 'omitnan');
std_dixot2  = std(iters_dixotomhsh2, 'omitnan');

fprintf("ΜΕΣΟΣ ΟΡΟΣ ΕΠΑΝΑΛΗΨΕΩΝ: %.2f\n", mean_dixot2);
fprintf("ΤΥΠΙΚΗ ΑΠΟΚΛΙΣΗ: %.2f\n", std_dixot2);


%διχοτόμηση klasikh
function [root, iterations] = bisection(f, a, b)
    tol = 1e-7;
    if f(a) * f(b) > 0
        warning('Δεν υπάρχει ρίζα στο [a,b]');
    end
    
    iterations = 0;
    while (b - a) / 2 > tol
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
    
    root = (a+b) / 2;  
end

%διχοτόμηση tropopoihmenh
function [root, iterations] = bisection2(f, a, b)
    tol = 1e-7;

    if f(a) * f(b) > 0
        error('Δεν υπάρχει ρίζα στο [a,b]');
    end

    iterations = 0;

    fa = f(a);
    fb = f(b);
    while (b-a)/2 > tol
        
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
            fb= fm;
        else         
            a = m;        
            fa = fm;
        end
    
        iterations = iterations+1;
    end
  
    root = (a+b) /2;
end