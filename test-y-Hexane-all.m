entA_x0=[];
entB_x0=[];
entR_x0=[];

difA_x0=[];
difB_x0=[];
difR_x0=[];

deA_x0=[];
deB_x0=[];
deR_x0=[];
current_dir = ['.\MyProjects\LIF-PIV\5MPa_Hexane_'];
for pamb=[40,70,100]
    for T=[30,60,90]
        
        for time = [3,4,5,6,7]
            figure
        %     time = 4;
            field=readimx([current_dir,num2str(pamb),'kPa_',num2str(T),'C\AvgV-flow\B',num2str(time,'%05d'),'.vc7']);
            [A,B,XData,YData]=showimx(field);
            close
            
            y=YData(27,:);
            x=XData(27,:);
        %     
        %     [y_max,y_max_index]=max(-y);
        %     [y_min,y_min_index]=min(-y);
            x=-y;
            [x_max,x_max_index]=max(-y);
            [x_min,x_min_index]=min(-y);
            [vor_x_val,vort]=min(abs(-y((6+3*time):x_max_index)));
            vort=vort+4+3*time;
            
            % X Fitting
            
            [dif_fit_result_x,dif_fit_eval_x]=diffusionFit(x(1:50),min(x_max_index,47));
            
            [ent_fit_result_x,ent_fit_eval_x]=entrainFit(x,max(7,x_min_index));
            
            vort=vort+1;
            
            
            ent_curve=[zeros(5,1);ent_fit_result_x(6:vort);zeros(51-(vort),1)]';
            dif_curve=[zeros(vort-1,1);flipud(dif_fit_result_x(x_max_index+1:x_max_index*2-vort+1));dif_fit_result_x(x_max_index+1:51)]';
            % dif_curve=[zeros(vort,1);dif_fit_result_x(vort+1:51)]';
            
            exp_curve=x-dif_curve-ent_curve;
            [de_fit_result_x,de_fit_eval_x]=deentrainFit(exp_curve,min(x_min_index,vort-3),vort);
        
        %    [de_fit_result_x,de_fit_eval_x]=deentrainFit(exp_curve,vort,x_min_index);
        %     
        %     x_min_index=0;
            exp_curve_r=[zeros(x_min_index-1,1);de_fit_result_x(x_min_index:vort-1);zeros(52-(vort),1)]';
            ent_curve=[zeros(4,1);ent_fit_result_x(5:vort-1);zeros(52-(vort),1)]';
            fitted_x=ent_curve+dif_curve+exp_curve_r(1:51);
            
            
            
            
            
            figure
            h5= plot(x,'-o','LineWidth',2);
            hold on
            h1= plot(fitted_x,'-','LineWidth',2);
            h2=plot(ent_curve,'--','LineWidth',1.2);
            h3=plot(dif_curve,'--','LineWidth',1.2);
            h4=plot(exp_curve_r,'--','LineWidth',1.2);
            
            legend('Experimental Data','Fitted Result','Entrainment','Diffusion','Recirculation')
            ylim([-10,10])
            
            
            % Fitting Result Extraction
            entA_x=ent_fit_result_x.a;
            entB_x=ent_fit_result_x.b;
            entR_x=ent_fit_eval_x.rsquare;
        
            difA_x=dif_fit_result_x.a;
            difB_x=dif_fit_result_x.b;
            difR_x=dif_fit_eval_x.rsquare;
        
            deA_x=de_fit_result_x.a;
            deB_x=de_fit_result_x.b;
            deR_x=de_fit_eval_x.rsquare;
            
            entA_x0=[entA_x0,entA_x]
            entB_x0=[entB_x0,entB_x];
            entR_x0=[entR_x0,entR_x];
            
            difA_x0=[difA_x0,difA_x];
            difB_x0=[difB_x0,difB_x];
            difR_x0=[difR_x0,difR_x];
            
            deA_x0=[deA_x0,deA_x];
            deB_x0=[deB_x0,deB_x];
            deR_x0=[deR_x0,deR_x];
            
            xlabel('y, mm')
            ylabel('Injection-induced Velocity $v_T$, $m/s$','interpreter','latex')
            
            title(['Fuel Temperature - ',num2str(T),'C ','0.',num2str(time),'ms  ','$P_{amb}$=',num2str(pamb),'KPa'],'interpreter','latex')
            % Fitting Result Presentation
        end
    end
end