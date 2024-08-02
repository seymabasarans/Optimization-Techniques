syms x1 x2 x3 x4
f = 10*(x2 - x1.^2)^2 + (1 - x1)^2 + 9*(x4 - x3.^2)^2 + (1 - x3)^2 + ((x2 - 1)^2 + (x4 - 1)^2) + 1.98*(x2 - 1)*(x4 - 1)
vars=[x1,x2,x3,x4];
hess = hessian(f,vars);
grad = gradient(f,vars);
% Error bound
epsilon = 10^-4;




%% Newton-Raphson
fprintf('Newton-Raphson Algorithm\n');

x0=-2+4*rand(4,1); %-->
%x=[-1;1]; %başlangıç vektörü
x=x0;
x0
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f\n',x(1),x(2),x(3),x(4),func(x))
x_next=x-inv(hessianfunc(x))* gradfunc(x);
fprintf('k=2, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),abs(func(x_next)-func(x)))
k=3;
%while(norm(gradfunc(x_next))>epsilon) normlu durdurma fonksiyonu
while(abs(func(x_next)-func(x))>epsilon)
   x=x_next;
   x_next=x-inv(hessianfunc(x))* gradfunc(x);
   fprintf('k=%d, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),abs(func(x_next)-func(x)))
   k=k+1;
end
toc
% title('Newton-Raphson Algorithm')
% set(gca,'fontsize',35)
% set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Hestenes-Stiefel Algorithm

fprintf('Hestenes-Stiefel Algorithm\n');
% x=rand(2,1);
% x=[-10;10];
x=x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f\n',x(1),x(2),x(3),x(4),func(x))
g=gradfunc(x);
d=-g;

% alpha argmin procedure
alpha=0:0.01:100000;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha);
alpha=alpha(ind);
% end of alpha argmin procedure
x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta=(g_next'*(g_next-g))/(d'*(g_next-g));
d_next=-g_next+beta*d;

% fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
fprintf('k=2, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),norm(gradfunc(x_next)))
k=3;
while(norm(gradfunc(x_next))>epsilon)
    x=x_next;
    g=g_next;
    d=d_next;

    % alpha argmin procedure
    alpha=0:0.01:100000;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    % end of alpha argmin procedure

    x_next=x+alpha*d;
    g_next=gradfunc(x_next);
    beta=(g_next'*(g_next-g))/(d'*(g_next-g));
    d_next=-g_next+beta*d;

    % fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    fprintf('k=%d, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',k,x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),norm(gradfunc(x_next)))
    k=k+1;
end
toc

%%  Polak-Ribiere  Algorithm

fprintf(' Polak-Ribiere  Algorithm\n');
% x=rand(2,1);
% x=[-10;10];
x=x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f\n',x(1),x(2),x(3),x(4),func(x))
g=gradfunc(x);
d=-g;

% alpha argmin procedure
alpha=0:0.01:100000;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha); 
alpha=alpha(ind);
% end of alpha argmin procedure
x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta = (g_next' * (g_next - g)) / (g' * g); 
d_next=-g_next+beta*d;

% fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
fprintf('k=2, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),norm(gradfunc(x_next)))
k=3;
while(norm(gradfunc(x_next))>epsilon)
    x=x_next;
    g=g_next;
    d=d_next;

    % alpha argmin procedure
    alpha=0:0.01:100000;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    % end of alpha argmin procedure

    x_next=x+alpha*d;
    g_next=gradfunc(x_next);
    beta = (g_next' * (g_next - g)) / (g' * g); 
    d_next=-g_next+beta*d;

    % fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    fprintf('k=%d, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',k,x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),norm(gradfunc(x_next)));
    k=k+1;
end
toc
%% Fletcher-Reeves Algorithm

fprintf('Fletcher-Reeves Algorithm\n');
% x=rand(2,1);
% x=[-10;10];
x=x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f\n',x(1),x(2),x(3),x(4),func(x))
g=gradfunc(x);
d=-g;

% alpha argmin procedure
alpha=0:0.01:100000;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha);
alpha=alpha(ind);
% end of alpha argmin procedure
x_next = x + alpha*d;
g_next = gradfunc(x_next);
beta = (g_next' * g_next) / (g' * g);  % Fletcher-Reeves
d_next = -g_next + beta * d;

% fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
fprintf('k=2, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),norm(gradfunc(x_next)))
k=3;
while(abs(func(x_next)-func(x))>epsilon)
    x=x_next;
    g=g_next;
    d=d_next;

    % alpha argmin procedure
    alpha=0:0.01:100000;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    % end of alpha argmin procedure

    x_next = x + alpha*d;
    g_next = gradfunc(x_next);
    beta = (g_next' * g_next) / (g' * g);  % Fletcher-Reeves
    d_next = -g_next + beta * d;

    % fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    fprintf('k=%d, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',k,x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),abs(func(x_next)-func(x)))
    k=k+1;
end
toc

%% Dai-Yuan Algorithm

fprintf('Dai-Yuan Algorithm \n');
% x=rand(2,1);
% x=[-10;10];
x=x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f\n',x(1),x(2),x(3),x(4),func(x))
g=gradfunc(x);
d=-g;

% alpha argmin procedure
alpha=0:0.01:100000;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha);
alpha=alpha(ind);
% end of alpha argmin procedure
x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta = (g_next' * g_next) / (d' * (g_next - g));  % Dai-Yuan
d_next=-g_next+beta*d;

% fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
fprintf('k=2, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),norm(gradfunc(x_next)))
k=3;
while(norm(gradfunc(x_next))>epsilon)
    x=x_next;
    g=g_next;
    d=d_next;

    % alpha argmin procedure
    alpha=0:0.01:100000;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    % end of alpha argmin procedure

    x_next = x + alpha*d;
    g_next = gradfunc(x_next);
    beta = (g_next' * g_next) / (d' * (g_next - g));  % Dai-Yuan
    d_next = -g_next + beta * d;

    % fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    fprintf('k=%d, x1=%f, x2=%f,x3=%f, x4=%f, f(x)=%f, error=%f\n',k,x_next(1),x_next(2),x_next(3),x_next(4),func(x_next),norm(gradfunc(x_next)))
    k=k+1;
end
toc