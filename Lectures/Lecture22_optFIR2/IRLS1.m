% m-file IRLS1.m to find the optimal solution to Ax=b % minimizing the L_p norm ||Ax-b||_p, using IRLS.
% Newton iterative update of solution, x, for M > N.
% For 2<p<infty, use homotopy parameter K = 1.01 to 2
% For 0<p<2, use K = approx 0.7 - 0.9
% csb 10/20/2012
function x = IRLS1(A,b,p,K,KK)
if nargin < 5, KK=10; end;
if nargin < 4, K = 1.5; end;
if nargin < 3, p = 10; end;
pk = 2;
x = pinv(A)*b;
E = [];
for k = 1:KK
	if p >= 2,
		pk = min([p, K*pk]);
	else
		pk = max([p, K*pk]);
	end
	e = A*x - b;
	w = abs(e).^((pk-2)/2);
	W = diag(w/sum(w));
	WA = W*A;
	x1 = (WA'*WA)\(WA'*W)*b;
	q = 1/(pk-1);
	if p>2,
		x=q*x1+(1-q)*x;
		nn=p;
	else
		x = x1;
		nn=2; 
	end
	ee = norm(e,nn);   
	E = [E ee];
end

figure,
plot(E)