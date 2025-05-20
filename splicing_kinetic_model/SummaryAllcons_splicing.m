[consgroups,consname_groups] = findgroups(sum_allwell(:,2));
fit_info=cell(size(consname_groups,1),8);
% fit_info :
% col1 : consname
% col2 3 4 5 : k_alpha k_s delta_s deltap_fixed
% col6 7 8 : R^2 RMSE reduced-chi-square
for i=1:size(consname_groups,1)
    sum_well_percons=sum_allwell(consgroups==i,:);
    NumberofTime=size(sum_well_percons,1);
    t=zeros(NumberofTime,1);
    P=zeros(NumberofTime,1);
    P_err=zeros(NumberofTime,1);
    S=zeros(NumberofTime,1);
    S_err=zeros(NumberofTime,1);

    for j=1:NumberofTime
        t(j)=str2double(sum_well_percons{j,3});
        P(j)=sum_well_percons{j,23}./10^(4); % reduced scale to E4
        P_err(j)=sum_well_percons{j,24}./10^(4);
        S(j)=sum_well_percons{j,26}./10^(4);
        S_err(j)=sum_well_percons{j,27}./10^(4);
    end
    % curve fitting
    deltap_fixed=sum_well_percons{1,21};
    mu_current=sum_well_percons{1,20};
    % if information is missing, skip this construct
    if deltap_fixed ==0 || mu_current ==0
        continue;
    end
    y=[P;S];
    y_err=[P_err;S_err];
    var=ones(length(y),1);
    lb = [];
    ub = [];
    fun = @(para,t)Splicing_model_fixeddecay(t,para,var,deltap_fixed);
    [x,resnorm] = lsqcurvefit(fun,[20 4 0.2],t,y./var,lb,ub);

    t_fit=[0:0.01:6.5]';
    var=ones(2*size(t_fit,1),1);
    y_fit=Splicing_model_fixeddecay(t_fit,x,var,deltap_fixed);
    P_fit=y_fit(1:length(t_fit),:);
    S_fit=y_fit(length(t_fit)+1:end,:);
    figure
    errorbar(t,P,P_err,"ro")
    hold on
    errorbar(t,S,S_err,"bo")
    hold on
    plot(t_fit,P_fit,'r-',t_fit,S_fit,'b-')
    xlabel('t/h')
    ylabel('Intensity (x10^{4} AU)')
    legend('[P]','[S]','[P]fitting','[S]fitting','Location', 'northwest')
    xlim([0 6.5])

    fit_info{i,1}=consname_groups{i,1};
    fit_info{i,2}=x(1);
    fit_info{i,3}=x(2);
    fit_info{i,4}=x(3);
    fit_info{i,5}=deltap_fixed;
    % goodness of fitting
    SS_res=resnorm;
    SS_total=sum(((y-mean(y)).^2));
    rsquare=1-SS_res./SS_total;
    RMSE=sqrt(resnorm./length(y));
    var=ones(2*size(t,1),1);
    y_pred=Splicing_model_fixeddecay(t,x,var,deltap_fixed);
    chi_square=sum(((y_pred-y)./(y_err)).^2);
    dof=length(y)-3;
    reduced_chi_square=chi_square./dof;

    fit_info{i,6}=rsquare;
    fit_info{i,7}=RMSE;
    fit_info{i,8}=reduced_chi_square;

    title(['construct = ' consname_groups{i,1}])
    f = gcf;
    exportgraphics(f,[summary_folder filesep consname_groups{i,1} '_mean.png'],'Resolution',300)
    close all



end

save([summary_folder filesep 'sum_allcons_mean.mat'],'fit_info')
