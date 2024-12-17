clear all;

rng(1);
rand_p = rand(3,100);
rand_p(1,:) = 1+rand_p(1,:);
rand_p(2,:) = 1+2*rand_p(2,:);
rand_p(3,:) = 1+2*rand_p(3,:);
e_dmd = zeros(1,100);
for i = 1:100
    e_dmd(i) = linear_stanford(rand_p(1,i),rand_p(2,i),rand_p(3,i));i
end

save('state100.mat',"e_dmd","rand_p");



