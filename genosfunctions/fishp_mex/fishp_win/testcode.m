clc;

%% Initialization
n_test = 5000;
mat_time = 0;
my_time = 0;

%% Main code
fprintf('Checking the correctness: ');
for ti=1:n_test
	org = [randn(5000,1)-rand()*3 randn(5000,1)+rand()*3 randn(5000,1)-rand()*3 randn(5000,1)+rand()*3]>0.5;
	X_hit = org(:,1) & org(:,2);
	Y_hit = org(:,3) & org(:,4);

	int_pxpy = sum( X_hit &  Y_hit);
	int_nxpy = sum(~X_hit &  Y_hit);
	int_pxny = sum( X_hit & ~Y_hit);
	int_nxny = sum(~X_hit & ~Y_hit);
	tbl = [int_pxpy int_nxpy; int_pxny int_nxny];
	
	tic;
	[~, mat_mute] = fishertest(tbl, 'Tail', 'left');
	[~, mat_coop] = fishertest(tbl, 'Tail', 'right');
	mat_time = mat_time + toc;
	
	tic;
	[my_mute, my_coop] = FastFisherExactTest(int_pxpy, int_nxpy, int_pxny, int_nxny);
	my_mute = exp(my_mute);
	my_coop = exp(my_coop);
	my_time = my_time + toc;
	
	if abs(mat_coop - my_coop)>mat_coop/1e10 || abs(mat_mute - my_mute)>mat_mute/1e10
		error('Results are not equal.'); 
	end 
end
fprintf('all tests went fine.\n');
fprintf('Mat time = %0.3fs\nMy time = %0.3fs\n', mat_time, my_time);
fprintf('Speed up: %0.1f times\n', mat_time/my_time);
