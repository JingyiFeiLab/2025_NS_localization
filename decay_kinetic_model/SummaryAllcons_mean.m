[consgroups,consname_groups] = findgroups(sum_allwell(:,2));
fit_info=cell(size(consname_groups,1),1);
% Column 1: consname  
%     - Name of the RNA construct being analyzed.
%
% Column 2: ln(I₀)  
%     - Fitted natural logarithm of the initial fluorescence intensity I₀.  
%       This is the y-intercept of the linearized decay model.
%
% Column 3: -δₚ  
%     - Fitted decay rate (delta_p), derived from the slope of the line fit to ln(I) vs. time.  
%       Since the exponential decay model is:  
%           I(t) = I₀ * exp(-δₚ * t),  
%       it is linearized as:  
%           ln(I(t)) = ln(I₀) - δₚ * t
%
% Column 4: R²  
%     - Coefficient of determination (R²), representing the goodness-of-fit of the linear model.  
%       Values closer to 1 indicate better fit.
for i=1:size(consname_groups,1)
    sum_well_percons=sum_allwell(consgroups==i,:);
    NumberofTime=size(sum_well_percons,1);
    t=zeros(NumberofTime,1);
    Pc=zeros(NumberofTime,1);
    Pc_se=zeros(NumberofTime,1);

    for j=1:NumberofTime
        t(j)=str2double(sum_well_percons{j,3});
        Pc(j)=sum_well_percons{j,10};
        Pc_se(j)=sum_well_percons{j,12};
    end

    lnPc=log(Pc);
    lnPc_se=Pc_se./Pc;

    [decay,gof_decay] = fit(t,lnPc,'poly1');
    coeffAll = coeffvalues(decay);

    fit_info{i,1}=consname_groups{i,1};
    fit_info{i,2}=coeffAll(1,2);
    fit_info{i,3}=-coeffAll(1,1);
    fit_info{i,4}=gof_decay.rsquare;

    t_fit=0:0.01:3.5;
    y_fit=coeffAll(1,1).*t_fit+coeffAll(1,2);
    figure
    plot(t_fit,y_fit,'b-')
    hold on
    errorbar(t,lnPc,lnPc_se,'o')
    xlabel('t /h')
    ylabel('ln[Pc]')
    title(['cons = ' consname_groups{i,1} ', ' '{\delta}_{p} = ' num2str(-coeffAll(1,1)) ', rsquare = ' num2str(gof_decay.rsquare)])
    f = gcf;
    exportgraphics(f,[summary_folder filesep consname_groups{i,1} '_mean.png'],'Resolution',300)
    close all



end

save([summary_folder filesep 'sum_allcons_mean.mat'],'fit_info')
