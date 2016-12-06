%%
clear;

%Step sizes to use for h
hVals = [0.1; 0.05; 0.01; 0.005; 0.001];

%Numerical error at each step size used
error = zeros(size(hVals));

%Loop through and test the step sizes
for i = 1:size(hVals)
    %Current step size
    h = hVals(i);
    %Time range
    time = (0 : h : 5);
    %Numerical integration values
    U = zeros(size(time))';
    %Half-step trapezoidal rule values
    Ustar = zeros(size(time))';
    %Initila value for the IVP
    U(1) = 1.5;
    t = 0;
    n = 1;
    lmda = -(10^6);

    %Run TR-BDF2 algorithm on given function
    while t < time(end)
        %Given function f(Un)
        fUn = lmda * (U(n) - cos(t)) - sin(t);
        %Intermediate trapezoidal step at t+h/2
        Ustar(n) = (1/(1-lmda*h/4)) * (U(n)...
            + (h/4)*(fUn - lmda * cos(t+h/2) - sin(t+h/2)));
        %BDF2 from t+h/2 --> t+h
        U(n+1) = (1/(3-lmda*h)) * (4*Ustar(n)...
            - U(n) - lmda*h*cos(t+h) - h*sin(t+h));
        
        n = n + 1;
        t = t + h;
    end
    
    fprintf('Numerically computed value for h = %.3f: u(t=5) = %.4f\n', h, U(end));

    %Global error with current step size
    error(i) = U(end) - 0.2837;
end

%Plot step size vs global error
plot(hVals, error);
snapnow