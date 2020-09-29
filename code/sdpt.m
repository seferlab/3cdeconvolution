cd('/Users/esefer/Downloads/SDPT3-4.0')
startup
[blk,At,C,b] = read_sdpa('/Users/esefer/Desktop/3cdeconvolution/code/tempsdpt');
[obj,X,y,Z] = sqlp(blk,At,C,b);
Xmat = cell2mat(X);
cd('/Users/esefer/Desktop/3cdeconvolution/code')
save('sdptoutput.mat','Xmat','obj')
