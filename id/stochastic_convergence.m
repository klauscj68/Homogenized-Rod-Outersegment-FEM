% Generates the figures of relative differences for reponses relevant to
% different living rhodopsins

clc
close all






load E_st_sal_III;


E=massa_E_st{1};

load E_st_sal_I;


E=[E;massa_E_st{1}];

load E_st_sal_II;


E=[E;massa_E_st{1}];

load E_st;


E=[E;massa_E_st{1}];


n_sample=size(E,1);

CV_area=zeros(1,n_sample);
CV_peak=zeros(1,n_sample);

E_area=trapz(time_downsample,E,2);
E_peak=max(E,[],2);


for cont=1:n_sample
    

    CV_area(cont)=std(E_area(1:cont))/mean(E_area(1:cont));
    CV_peak(cont)=std(E_peak(1:cont))/mean(E_peak(1:cont));

end
figure(2)
plot(1:n_sample,CV_area,'b-');
a=title('CV of the E^* area vs the sample number. Simulations up to 8s');
set(a,'fontsize',18);
set(gca,'fontsize',18);
a=xlabel('sample number');
b=ylabel('CV of the E^* area');
set([a b],'fontsize',20);
figure(1)
plot(1:n_sample,CV_peak,'b-');
a=title('CV of the peak of E^* vs the sample number. Simulations up to 8s');
set(a,'fontsize',18);
set(gca,'fontsize',18);
a=xlabel('sample number');
b=ylabel('CV of the E^* peak');
set([a b],'fontsize',20);

