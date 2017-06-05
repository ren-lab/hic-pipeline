clear
clc

dir=getenv('TAD_PATH')
addpath(genpath(strcat(dir,'/required_modules/HMMall/'))); 
addpath(genpath(strcat(dir,'/required_modules/pmtk3-2april2011/'))); 
addpath(genpath(strcat(dir,'/required_modules/order')));
input = getenv('TAD_INPUT')
output = getenv('TAD_OUTPUT')

%addpath(genpath('required_modules/HMMall/')); 
%addpath(genpath('required_modules/pmtk3-2april2011/')); 
%addpath(genpath('required_modules/order'));

minaic = 100000000;
filename = input;
%data = load(filename);
data=dlmread(filename,'\t',0,3);
prior0    = [0.33 0.33 0.33];
transmat0 = [0.33 0.33 0.33
             0.33 0.33 0.33
             0.33 0.33 0.33];
testX = data(:,1)';
endm  = 20;
aic   = ones(1,endm)*10000000;    
inim  = 1;

for M=inim:endm
    M
       
    Q = 3; O = 1;
    [mixModel] = mixGaussFit(testX',Q*M,'maxIter', 500);
    mu0     = mixModel.cpd.mu;
    Sigma0  = bsxfun(@plus, eye(1), mixModel.cpd.Sigma);
    mixmat0 = mixModel.mixWeight;
    mu0     = reshape(mu0, [O Q M]);
    Sigma0  = reshape(Sigma0, [O O Q M]);
    mixmat0 = reshape(mixmat0, [Q,M]);
        
    [LL, prior1, transmat1, mu1, Sigma1, mixmat1]     = mhmm_em(testX, prior0, transmat0, mu0, Sigma0, mixmat0,'max_iter', 500,'thresh',1e-5);

    loglik = mhmm_logprob(testX, prior1, transmat1, mu1, Sigma1, mixmat1);
    B = mixgauss_prob(testX, mu1, Sigma1, mixmat1);
    path = viterbi_path(prior1, transmat1, B);
    [gamma, alpha, beta, logp] = hmmFwdBack(prior1, transmat1, B);
    decodedFromEMmaxMarg = maxidx(gamma);

    s_p = size(prior1);s_t = size(transmat1);s_m = size(mu1);s_v = size(Sigma1);s_i = size(mixmat1);s_te = size(testX);
        
    if (M==1)
        num_p = (s_p(1)*s_p(2)) + (s_t(1)*s_t(2)) + (s_m(1)*s_m(2)) + (s_v(1)*s_v(2)) + (s_i(1)*s_i(2));
    else
        num_p = (s_p(1)*s_p(2)) + (s_t(1)*s_t(2)) + (s_m(1)*s_m(2)*s_m(3)) + (s_v(1)*s_v(2)*s_v(3)*s_v(4)) + (s_i(1)*s_i(2));
    end     
        
    aic(M)    = (-2*loglik) + (2*num_p);
    ll(M)     = loglik;
    best_p(M,:) = path';
    best_d(M,:) = decodedFromEMmaxMarg';
    best_g1(M,:) = gamma(1,:)';
    best_g2(M,:) = gamma(2,:)';
    best_g3(M,:) = gamma(3,:)';
end
          
ord = order(min(aic)) - 1;
div = power(10,ord);

for k=inim:endm
     min(aic);
     paic(k) = exp((min(aic)-aic(k))/(div*2));
     if ((paic(k)) >= 0.9)
         ind = k;
         break;
     end
end    
        
final_p      = best_p(ind,:)';
final_d      = best_d(ind,:)';
final_g(:,1) = best_g1(ind,:)';
final_g(:,2) = best_g2(ind,:)';
final_g(:,3) = best_g3(ind,:)';
final_aic = aic(ind); minaic = min(aic);

states = [final_p final_d final_g];
file =  output;
save (file,'aic','ind','states','-ASCII');      

exit;
