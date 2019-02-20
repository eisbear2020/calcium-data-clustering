function mi = mutualInfoHist(M)

   [nR, nC] = size(M);

   eps = 1e-7;

   mi = 0;

   jointP = M ./ sum(M(:));

   probRow = sum(jointP, 2);
   probCol = sum(jointP, 1);

   for i = 1:nR
       iP = probRow(i);
       for j = 1:nC
           jP = probCol(j);
           if jointP(i,j) > eps
               mi = mi + (jointP(i,j) * (log2(jointP(i,j)) -  log2(iP) - log2(jP)));
           end
       end
   end
   
end
