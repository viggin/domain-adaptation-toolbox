 function plotpra(PRA,index)
 

 FIT=PRA.FIT;
 PRED=PRA.PRED;
 [M,N]=size(FIT);
 rsf=FIT(:,index);rsf=rsf(find(isnan(rsf)==0));
 rsp=PRED(:,index);rsp=rsp(find(isnan(rsp)==0));
 FIT=reshape(FIT,M*N,1);
 PRED=reshape(PRED,M*N,1); 
 FIT=FIT(find(isnan(FIT)==0));
 PRED=PRED(find(isnan(PRED)==0));
 
 % plot
 nbin=30;
 subplot(221)
 hist(FIT,nbin);
 subplot(222)
 hist(PRED,nbin);
 subplot(223)
 hist(rsf,nbin);
 subplot(224)
 hist(rsp,nbin);
 