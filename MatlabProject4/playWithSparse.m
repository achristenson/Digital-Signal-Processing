% this shows an example of using sparse matrices in matlab to save space

% let's define a smallish (in case you want to look at it)
% matrix of random numbers

mtxRandom = rand(5,4) 
% leaving the semicolon off because it's so much more fun that way!! 

disp('pausing; hit any key to continue');
pause

% now, set any values below a threshold to zero
mtxZeroed = mtxRandom;
mtxZeroed(mtxZeroed<0.1) = 0

disp('pausing; hit any key to continue');
pause
% now use the 'sparse' command to compress the matrix
mtxSparse = sparse(mtxZeroed)