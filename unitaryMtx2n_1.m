function Q = unitaryMtx2n_1(n)
%Computes a left
%INPUT:
%   n: size of output matrix will be (2*n + 1) x (2*n + 1)
%OUTPUT:
%   Left real unitary matrix of size 2n x 2n
Q = (1/sqrt(2))*[eye(n) zeros(n,1) j*eye(n);
      zeros(n,1)' sqrt(1/2) zeros(n,1)';
        flip(eye(n)) zeros(n,1) -j*flip(eye(n))];
        
%remove center row and column
% [row, column] = size(Q);
% Q(round(row/2),:) = [];
% Q(:,round(column/2)) = [];

end