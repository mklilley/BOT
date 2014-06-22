%     Copyright 2011 Matthew Lilley

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License Version 3 
%     as published by the Free Software Foundation 
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>


%    BOT v1.9 06-01-2011 Created by Matthew Lilley.  

%    It is the author's wish to be included as an author on the first two 
%    publications written using results obtained from this code.

%    The author would be grateful for any suggestions regarding improvement of the code.
%    Please contact Matthew Lilley at bumpontail@gmail.com.  


%This code simulates a single sinusoidal mode in the bump on tail problem.
%The 'spatial' grid here actually has units of time.  It is like
%1/(kv-oemga).  Frequency is normalised to gamma_L and and time is
%normalised to 1/gamma_L

%To compare results of this code to Berk-Breizman PRL one should convert
%units of time from gamma_L to units of gamma_L-gamma_d by multiplying by 
%(1-gamma_d/gamma_L).  One should also convert units of bounce frequency 
%squared to units of A by dividing by (1-gamma_d/gamma_L)^(5/2)

clear all

%Determine if your computer is a PC or not, for visualisation output
PC=ispc;

fprintf(1, 'BOT v1.9.1 13-01-2014 Created by Matthew Lilley\n')

%obtains simulations parameters from user input

descrip_value={'number of conjugate velocity points (MUST be an odd number)','size of computational box in units of 1/\gamma_L','time step/spatial step (MUST be an integer)','Simulation time in units of 1/\gamma_L','Maximum harmonic number F e.g 3 gives f_0, f_1, f_2 and f_3 (MIN=1)','Number of loop iterations',...
    'Initial bounce frequency in units of \gamma_L','\gamma_d in units of \gamma_L','Krook collision frequency in units of (\gamma_L-\gamma_d)','Diffusive collision frequency in units of (\gamma_L-\gamma_d)','Drag collision frequency in units of (\gamma_L-\gamma_d)','Visualisation (0=off, 1=default, 2=change settings)',...
    'Timing options (1=all effects on continuously, 2=change options)'};
    


f = figure('WindowStyle','modal','Visible','off','units','normalized','Position',[0.1,0.1,0.8,0.8],'menubar','none');
AxesHandle=axes('Parent',f,'Position',[0,0,1,1],'Visible','off');

for j=1:8
user(j) = uicontrol('Style','edit','units','normalized','Position',[0.25-0.125,0.9-(j-1)*0.1,0.25,0.05]);
descrip(j)=text(0.25-0.125,0.975-(j-1)*0.1,descrip_value(j),'Interpreter','tex','units','normalized','Parent',AxesHandle);
end

for j=9:13
user(j) = uicontrol('Style','edit','units','normalized','Position',[0.75-0.125,0.9-(j-9)*0.1,0.25,0.05]);
descrip(j)=text(0.75-0.125,0.975-(j-9)*0.1,descrip_value(j),'Interpreter','tex','units','normalized','Parent',AxesHandle);

end

h = uicontrol('units','normalized','Position',[0.5-0.1,0.1,0.2,0.1],'String','Continue',...
              'Callback','uiresume(gcbf)');

default={'1001','10','1','300','10','2','1e-3','0.9','0','0.3','0','1','1'};
    



for j=1:13
set(user(j),'String',default(j));
end
set(f,'Visible','on');

uiwait(gcf); 


answer=get(user,'string');

close(f);
clear f user descrip h AxesHandle

if str2num(answer{13}{1})==1
    answer{13}{1}=answer{4}{1};
    answer{14}{1}=answer{4}{1};
    answer{15}{1}=answer{4}{1};
    answer{16}{1}=answer{4}{1};
    answer{17}{1}=answer{4}{1};
else
    
    descrip_value={'Equilibrium slope on until t=...','\gamma_d on until t=...','\beta on until t=...','\nu on until t=...','\alpha on until t=...'};
    f = figure('WindowStyle','modal','Visible','off','units','normalized','Position',[0.425,0.1,0.15,0.8],'menubar','none');
AxesHandle=axes('Parent',f,'Position',[0,0,1,1],'Visible','off');

for j=1:5
user(j) = uicontrol('Style','edit','units','normalized','Position',[0.5-0.45,0.7-(j-1)*0.1,0.9,0.05]);
descrip(j)=text(0.5-0.45,0.775-(j-1)*0.1,descrip_value(j),'Interpreter','tex','units','normalized','Parent',AxesHandle);
end



h = uicontrol('units','normalized','Position',[0.5-0.45,0.1,0.9,0.1],'String','Continue',...
              'Callback','uiresume(gcbf)');
          
          default={answer{4}{1},answer{4}{1},answer{4}{1},answer{4}{1},answer{4}{1}};
          
          for j=1:5
set(user(j),'String',default(j));
end
set(f,'Visible','on');

uiwait(gcf); 


answer2=get(user,'string');

    answer{13}{1}=answer2{1}{1};
    answer{14}{1}=answer2{2}{1};
    answer{15}{1}=answer2{3}{1};
    answer{16}{1}=answer2{4}{1};
    answer{17}{1}=answer2{5}{1};

close(f);
clear f user descrip h AxesHandle
    
end
    
 

%Starts the clock to time the simulation
tic

%Records the start time as a stamp for any saved files
cl=clock;

%number of 'spatial points' - This must be an odd number to put the zero in
%the centre of the grid (1001 is a good number to start with)
npx= str2num(answer{1}{1});
%Represent npx variable as a string that can be printed in graphs etc
npxs=num2str(npx);

%defines the 'spatial grid' (10 is a good number to start with for pbox)
pbox= str2num(answer{2}{1});
p=linspace(-pbox,pbox,npx);
%Represent pbox variable as a string that can be printed in graphs etc
pboxs=num2str(pbox);

% determines 'spatial' resolution
dp=p(2)-p(1);

%calculates the position of the zero
p0=(npx-1)/2+1;

%Corrects something wrong with MATLAB linspace routine.  The issue is that
%using linspace doesnt give you an antisymmetric p vector, there are small
%errors. These propagate and mess up the IFFT, basically you end up with an
%complex f_0 if you dont correct for this.
p(1:p0-1)=-fliplr(p(p0+1:end));

%Define the temporal resolution - the time step must be an integer multiple
%of the 'spatial step', this means advection is performed perfectly with no
%interpolation
nt= str2num(answer{3}{1});
%Represent nt variable as a string that can be printed in graphs etc
nts=num2str(nt);
dt=dp*nt;
%Represent dt variable as a string that can be printed in graphs etc
dts=num2str(dt);
Freq=1/dt;
%Defines the 'frequency in p space', this is used when converting F back to
%real space
Freqp=1/dp;

%Defines simulation time
tend= str2num(answer{4}{1});

%Determines what kind of visualisation the user selected
question=str2num(answer{12}{1});

%calculates the total number of time steps rounded to the nearest fram_num (the
%simulation is split into fram_num pieces either chosen by the user or 100 as default)
if question==2
    frame_num=str2num(cell2mat(inputdlg('How many frames do you want in the visualisation?','Frames',1,{'100'})));
else
    frame_num=100;
end
maxi=frame_num*round(tend/dt/frame_num);
%creates the time vector
time=[0:maxi]*dt;

%Defines the total number of harmonics. The minimum harmonic number
%must be at least 1
while str2num(answer{5}{1}) < 1
    fprintf(1,'\nERROR: Minimum harmonic number is 1 \n');
    nx=input('nx=...   ')+1;
end
nx= str2num(answer{5}{1})+1;

%Represent nx variable as a string that can be printed in graphs etc
nxs=num2str(nx-1);

%Defines the number of iterations used in the calculation of the distribution
%function: 2 seems to be enough.
loopi= str2num(answer{6}{1});

%Creates and empty 'electric field' vector and initialises the 'electric
%field' (1e-6 is good to start with), note that this is actually the square
%of the bounce frequency in units of gamma_L^2
E=zeros(1,maxi+1);
E(1)= str2num(answer{7}{1})^2;


%Creates empty vectors which are used in the calcualtion of the perturbed
%distribution function in fourier space, i.e 1/(kv-omega) space.
F=zeros(nx,size(p,2));
F2=F;
G=F;
F3=F;
F4=F;


%Creates a left shifted p vector, which is required in the advection process.
%The shift is greater for the higher harmonics.  Nothing is advected from
%outside the p box (from the right), so the shifted vector is padded with 
%zeros at the right hand side.
for llp=1:nx
pshift(llp,:)=[p(1+nt*(llp-1):end),zeros(1,nt*(llp-1))];
end


%Creates a special vector to deal with the contribution to the distrbution
%function from the equilibium part.  This contribution comes in the form of 
%a delta function at p=0, i.e at the point p(p0).
p2=zeros(1,npx);
p20=zeros(1,npx);
if nt==1
    p20(p0)=0.5;
    p20(p0-nt)=0.5;
else
    p20(p0)=0.5;
    p20(p0-nt+1:p0-1)=1;
    p20(p0-nt)=0.5;
end
p20=p20/sqrt(2*pi);


%Define the damping in untis of gamma_L
gd= str2num(answer{8}{1});
%Represent damping variable as a string that can be printed in graphs etc
gds=num2str(gd);

%Define the Krook, drag and diffusion collision frequencies in untis of 
%(gamma_L-gamma_d) 
betai= str2num(answer{9}{1});
nui= str2num(answer{10}{1});
alphai= str2num(answer{11}{1});

%Represent collision frequencies as strings that can be printed in graphs etc
betais=num2str(betai);
nuis=num2str(nui);
alphais=num2str(alphai);

%Calculate the Krook drag and diffusion collision frequencies in untis of 
%gamma_L
beta=betai*(1-gd);
nu=nui*(1-gd);
alpha=alphai*(1-gd);

%Calculate the collisional and electric field exponential factors required
%to produce the distribution function and electric field and the at the 
%next timestep.

%This exg vector turns off the damping at a time t=damptime
damptime=str2num(answer{14}{1});
exg=[exp(-gd*dt).*ones(1,floor(damptime/dt)),ones(1,size(time,2)-floor(damptime/dt))];
%note if you want the damping to turn off slowly you should use something
%like:
%exg=[exp(-gd*dt).*ones(1,floor(damptime/dt)),linspace(exp(-gd*dt),exp(-1.1*dt),4000),exp(-1.1*dt)*ones(1,size(time,2)-4000-floor(damptime/dt))];
%exg=[exp(-gd*dt).*ones(1,floor(damptime/dt)),linspace(exp(-gd*dt),exp(-1.1*dt),size(time,2)-floor(damptime/dt))];

exb=exp(-beta*dt);
for k=1:nx
    exn(k,:)=exp(-nu^3*dt/3*(3*p.^2+(3*p*(k-1)*dt)+((k-1)*dt)^2));
end
for k=1:nx
    exa(k,:)=exp(alpha^2*dt*i/2*(2*p+((k-1)*dt)));
end
EXPALL=exb.*exn.*exa;

%Defines the time at which Krook, diffusion and drag will turn off
krooktime=str2num(answer{15}{1});
difftime=str2num(answer{16}{1});
dragtime=str2num(answer{17}{1});

%Creates an exponential factor which forms part of the contribution from
%the delta function arising due to equilibrium part of the distribution
%function when collisions are present.
exbdel=zeros(1,npx);
for lll=p0-nt:p0
    exbdel(lll)=exp(beta*p(lll))*exp(nu^3/3*p(lll)^3)*exp(-i*alpha^2/2*p(lll)^2);
end



%Creates an array which turns the original slope of F_0 off at a time t=restime
restime=str2num(answer{13}{1});
resonoff=[ones(1,floor(restime/dt)),zeros(1,size(time,2)-floor(restime/dt))];



% If you choose to have non default visualisation then you get asked some
% questions

if question==1
    question0='Yes';
    question2='Yes';
    question3='Yes';
    question4='Real space   ';
elseif question==2
    options.Interpreter = 'tex';
    options.Default = 'Yes';

    question0=questdlg('Do you want to see a Fourier spectrogram, n.b this will appear at the end of the simulation','Spectrogram','Yes','No ',options);
    questionCont=questdlg('Do you want to see a phase space plot, n.b if you want this you will not be able to view anything else in real time','Phase space plot','Yes','No ',options);
    
    if questionCont=='Yes'
        %Defines things that are needed to make the phase space plot
        
        %Defines the upper and lower velocity limits, in units of
        %gamma_L/2/pi. 
         options.Interpreter = 'tex';
         options.Default = 'Yes';
         
        
        DEF={'10','-10','4*pi'};    
        prompt={'Upper "velocity" limit in units of \gamma_L','Lower "velocity" limit in units of \gamma_L','Spatial range'};
        A=inputdlg(prompt,'Phase space options',1,DEF,options);
        
        fv=str2num(A{1});
        sv=str2num(A{2});
        
        %Sets the minimum and maximum values that the contour plots of the
        %phase space will show
  
            Cmin=str2num(A{2})/pi;
            Cmax=str2num(A{1})/pi;
        
        x_max=str2num(A{3});
            
        %Tries to load a high resolution colour map if the user has it on
        %the system.  If not the code wont crash, it will just be ignored
        try load('mapn')
        catch exception
        end
                       
        %Creats a high resolution velocity vector    
        nU=8*npx+npx;
        O=linspace(-Freqp/2,Freqp/2,nU)*2*pi;
        start=find(abs(O-sv)==min(abs(O-sv)));
        finish=find(abs(O-fv)==min(abs(O-fv)));
        O=O(start:finish);
        
         x=[0:0.01:x_max];
        [xx,OO]=meshgrid(x,O);
        
        
        
         question2='No ';
         question3='No ';
         
    else
    
    question2=questdlg('Do you want to see the evolution of the bounce frequency','Bounce Frequency','Yes','No ',options);
    
    question3=questdlg('Do you want to see the evolution of F','Distribution Function','Yes','No ',options);
    if question3=='Yes'
        options.Default = 'Fourier space';
        question4=questdlg('How do you want to see F','Distribution Function','Fourier space','Real space   ',options);
        if question4=='Fourier space'
            question5=str2num(cell2mat(inputdlg('Which fourier component do you want to see (0,1,2,3,4,5 etc)','Distribution Function',1,{'1'})));
            if question5>nx-1
                question5=nx-1;
            end
            question5s=sprintf('%1.0f ',question5);
            options.Default = 'Real     ';
            question6=questdlg('Which part of F do you want to see','Distribution Function','Real     ','Imaginary',options);
        end
    end

    end
end
    
if question==1 | question==2
    %Initialises frames for movie production
    axis tight
    set(gca,'nextplot','replace');
    %This finds out whether you have a good compressor to create nice avi
    %files.  If you dont then it creates an uncompressed file.  In UNIX
    %only uncompressed avi files can be made. Type mmcompinfo into matlab
    %command prompt to look at what compressors you have.

    if PC==1
        qmpeg=mmcompinfo('video','MPEG-4');
        qcinepak=mmcompinfo('video','Cinepak');
        if qmpeg==-1 & qcinepak==-1
            aviobj=avifile(['BOT_a=' alphais '_n=' nuis '_b=' betais '_gd=' gds '_pmax=' num2str(loopi) '_N=' num2str(nx-1) '_smax=' num2str(pbox) '_spoints=' num2str(npx) '_(' num2str(cl(1)) '_' num2str(cl(2)) '_' num2str(cl(3)) '_' num2str(cl(4)) '_' num2str(cl(5)) ').avi'],'compression','NONE','fps',5,'quality',100)   
        elseif qmpeg==-1
            aviobj=avifile(['BOT_a=' alphais '_n=' nuis '_b=' betais '_gd=' gds '_pmax=' num2str(loopi) '_N=' num2str(nx-1) '_smax=' num2str(pbox) '_spoints=' num2str(npx) '_(' num2str(cl(1)) '_' num2str(cl(2)) '_' num2str(cl(3)) '_' num2str(cl(4)) '_' num2str(cl(5)) ').avi'],'compression','Cinepak','fps',5,'quality',100)
        else
            aviobj=avifile(['BOT_a=' alphais '_n=' nuis '_b=' betais '_gd=' gds '_pmax=' num2str(loopi) '_N=' num2str(nx-1) '_smax=' num2str(pbox) '_spoints=' num2str(npx) '_(' num2str(cl(1)) '_' num2str(cl(2)) '_' num2str(cl(3)) '_' num2str(cl(4)) '_' num2str(cl(5)) ').avi'],'compression','MPG4','fps',5,'quality',100)
        end
        
    elseif PC==0
        aviobj=avifile(['BOT_a=' alphais '_n=' nuis '_b=' betais '_gd=' gds '_pmax=' num2str(loopi) '_N=' num2str(nx-1) '_smax=' num2str(pbox) '_spoints=' num2str(npx) '_(' num2str(cl(1)) '_' num2str(cl(2)) '_' num2str(cl(3)) '_' num2str(cl(4)) '_' num2str(cl(5)) ').avi'],'compression','None','fps',5,'quality',100)
    end

end



%This loop splits the computation into 100 pieces, so that a movie frame can
%be produced at regular intervals, but not after every timestep which would
%slows down the code DRAMATICALLY.

gr=maxi/frame_num;
for gr2=0:(frame_num-1)
    gr2+1
    
    for j=1+gr2*gr:(gr2+1)*gr
       
        %This allows collisions to be turned off at a time chosed by the
        %user
        if j>=floor(krooktime/dt)
            exb=1;
            beta=0;
        end
        if j>=floor(difftime/dt)
            exn=ones(nx,npx);
            nu=0;
        end
        if j>=floor(dragtime/dt)
            exa=ones(nx,npx);
            alpha=0;
        end
                   
        EXPALL=exb.*exn.*exa;
            for lll=p0-nt:p0
            exbdel(lll)=exp(beta*p(lll))*exp(nu^3/3*p(lll)^3)*exp(-i*alpha^2/2*p(lll)^2);
            end
        
            
        %Stores F at the last time step
        G=F;
        
        
        %Advects F at the pervious time step to the left, the advection is
        %proportional to the harmonic number.  This will form part of the
        %first guess at the solution for F at the next time step
        for k=2:size(F,1)
            F(k,:)=[G(k,1+nt*(k-1):end),zeros(1,nt*(k-1))];
        end
        
        %The advection process also requires the collisional exponential
        %factors
        F=F.*EXPALL;
        
        %Calculates the first guess of the electric field at the next time
        %step by neglecting the non-linear terms
        E2=(0.5*dt*2*sqrt(2*pi)*(F(2,p0)+G(2,p0)*exg(j))+E(j)*exg(j))/(1-resonoff(j)*0.5*dt);
   
               
        %Calculates the first guess of F at the next time step, which is
        %just the free streaming solution except for F_1, which needs a
        %contribution from the delta function
        p2=p20.*[zeros(1,p0-nt-1),linspace(E(j),E2,nt+1),zeros(1,(npx-1)/2)];
        %You have to do this next line because of something strange
        %with the accuracy of linspace
        p2(p0)=0.5*E2/sqrt(2*pi);  
        F3=F;        
        F3(2,:)=F(2,:)+p2.*exbdel*resonoff(j);

        
        
        %Calculates the part of the non-linearity which only depends on
        %information from the previsous time step.  For F_1 for example you need to 
        %include contribution from F_2 and F_0 for a single harmonic of the 
        %electric field.  Again there is padding with zeroes when the advection 
        %runs out of known data points.
        
        %For F_0 the contribution is from F_1 and F_-1.  F_-1 is related to
        %F_1 by F_-1 = complexconjugate F_1(-p) hence the matrix flip.
        k=1;
        F2(k,:)=(0.5*dt*(i/2)*p.*(E(j)*fliplr(conj(G(k+1,:)))+conj(E(j))*G(k+1,:))).*EXPALL(k,:);
        
        for k=2:size(F,1)-1
            F2(k,:)=(0.5*dt*(i/2)*pshift(k,:).*(E(j)*([G(k-1,1+nt*(k-1):end),zeros(1,nt*(k-1))])+conj(E(j))*([G(k+1,1+nt*(k-1):end),zeros(1,nt*(k-1))]))).*EXPALL(k,:);
        end
        
        %The last fourier harmonic can only receive contributions from the
        %harmonic below it, as there is no higher harmonic left in the code
        k=size(F,1);
        F2(k,:)=(0.5*dt*(i/2)*pshift(k,:).*(E(j)*([G(k-1,1+nt*(k-1):end),zeros(1,nt*(k-1))]))).*EXPALL(k,:);
 
        
        %This iterative loop is to determine F and E at the next time step, which
        %we do not know.  We begin with a guess which is E2 and F3
            
        
        for l=1:loopi
            
                 
        %Calculates the part of the non-linearity which depends on
        %information from the NEXT time step, which needs to be iteratively
        %obtained
            k=1;
            F4(k,:)=0.5*dt*(i/2)*p.*(E2*fliplr(conj(F3(k+1,:)))+conj(E2)*F3(k+1,:));
            
            
            for k=2:size(F,1)-1
                F4(k,:)=0.5*dt*(i/2)*p.*(E2*(F3(k-1,:))+conj(E2)*F3(k+1,:));
            end
            
            k=size(F,1);
            F4(k,:)=0.5*dt*(i/2)*p.*(E2*(F3(k-1,:)));

          %Calculates the total non-linear contribution  
          F4=F4+F2;
          
          
            %Calculates the next estiamte for E and F at the next timestep
            %including the non-lineairty
            E2=(0.5*dt*2*sqrt(2*pi)*(F(2,p0)+F4(2,p0)+G(2,p0)*exg(j))+E(j)*exg(j))/(1-resonoff(j)*0.5*dt);
            %The contribtuion from the delta function must also be updated
            %whch the updated estiamte for E2
            p2(p0-nt:p0)=p20(p0-nt:p0).*linspace(E(j),E2,nt+1);
            %You have to do this next line because of something strange
            %with the accuracy of linspace
            p2(p0)=0.5*E2/sqrt(2*pi);            
            % The next guess for F is the free streaming function + the
            % non-linear correction + the delta functin piece.
            F3=F4+F;
            F3(2,:)=F3(2,:)+p2.*exbdel*resonoff(j);

            
 
        end
        
        %updates F and E ready for the next timestep
        F=F3;
        E(j+1)=E2;
        

%Checks to see if we have lost numerical stability in the iterative loop
if p(end)*dt*abs(E2)/2>1
    fprintf(1,'\nLOST STABILITY\n');
    break
end
             
        
    end

    

%Creates a label for the time that can be used in movie plots
timel=j*dt;
timels=sprintf('%1.1f ',timel);
    

%This section performs all the plots based on the visualiastion decisions
%made by the user.
if question==1 | question==2
    
    %Multiple plots are made in this case
    if question2=='Yes' & question3=='Yes'
        
        subplot(2,1,1)
        
        %Plots the Fourier transorm of F
        if question4=='Fourier space'
            if question6=='Real     '
                
                plot(p,real(F(question5+1,:)))
                
                title(['\beta/(\gamma_L-\gamma_d)=' betais ' ' '\nu/(\gamma_L-\gamma_d)=' nuis '  ' '\alpha/(\gamma_L-\gamma_d)=' alphais '  ' '\gamma_d/\gamma_L=' gds '  ' 't\times\gamma_L=' timels])
                xlabel('s (Fourier conjugate to (kv-\omega)/\gamma_L)')
                ylabel(['Re(f_{' question5s '})'])
            
            
            elseif question6=='Imaginary'
                
                plot(p,imag(F(question5+1,:)))
                
                title(['\beta/(\gamma_L-\gamma_d)=' betais ' ' '\nu/(\gamma_L-\gamma_d)=' nuis '  ' '\alpha/(\gamma_L-\gamma_d)=' alphais '  ' '\gamma_d/\gamma_L=' gds '  ' 't\times\gamma_L=' timels])
                xlabel('s (Fourier conjugate to (kv-\omega)/\gamma_L)')
                ylabel(['Im(f_{' question5s '})'])
            
            
            end
        
        %Plots f_0 in real space
        elseif question4=='Real space   '
            
            
            %Plots the spatially averaged part of distribution.
            plot(linspace(-Freqp/2,Freqp/2,10*npx+npx)*2*pi,real((fftshift(ifft(ifftshift([zeros(1,5*npx),F(1,:),zeros(1,5*npx)]))))*(10*npx+npx)*dp/sqrt(2*pi))+1/pi*2*pi*linspace(-Freqp/2,Freqp/2,10*npx+npx))
            % To plot it with without the equilibrium you should use:
            %plot(linspace(-Freqp/2,Freqp/2,10*npx+npx)*2*pi,real((fftshift(ifft(ifftshift([zeros(1,5*npx),F(1,:),zeros(1,5*npx)]))))*(10*npx+npx)*dp/sqrt(2*pi)))
            xlim([-10,10])
            ylim([-10/pi,10/pi])
            
            title(['\beta/(\gamma_L-\gamma_d)=' betais ' ' '\nu/(\gamma_L-\gamma_d)=' nuis '  ' '\alpha/(\gamma_L-\gamma_d)=' alphais '  ' '\gamma_d/\gamma_L=' gds '  ' 't\times\gamma_L=' timels])
            xlabel('(kv-\omega)/\gamma_L')
            ylabel('F_0 + f_0')
            
            
            
        end
        
        %Plots the bounce frequency
        subplot(2,1,2)
        plot(time(1:j),abs(E(1:j)))
        
        xlabel('t\times \gamma_L')
        ylabel('|\omega_B^2|/\gamma_L^2')
        
        
        FRAMEDAT=getframe(gcf);
        aviobj = addframe(aviobj,FRAMEDAT);
        
    %Only bounce frequency is plotted here    
    elseif question2=='Yes' & question3=='No '
        
        plot(time(1:j),abs(E(1:j)))
        title(['\beta/(\gamma_L-\gamma_d)=' betais ' ' '\nu/(\gamma_L-\gamma_d)=' nuis '  ' '\alpha/(\gamma_L-\gamma_d)=' alphais '  ' '\gamma_d/\gamma_L=' gds '  ' 't\times\gamma_L=' timels])
        xlabel('t\times \gamma_L')
        ylabel('|\omega_B^2|/\gamma_L^2')
        
        FRAMEDAT=getframe(gcf);
        aviobj = addframe(aviobj,FRAMEDAT);
    
    %Only the distribution function is plotted here    
    elseif question2=='No ' & question3=='Yes'
         
        %Plots the Fourier transorm of F
        if question4=='Fourier space'
            if question6=='Real     '
                
                plot(p,real(F(question5+1,:)))
                
                title(['\beta/(\gamma_L-\gamma_d)=' betais ' ' '\nu/(\gamma_L-\gamma_d)=' nuis '  ' '\alpha/(\gamma_L-\gamma_d)=' alphais '  ' '\gamma_d/\gamma_L=' gds '  ' 't\times\gamma_L=' timels])
                xlabel('s (Fourier conjugate to (kv-\omega)/\gamma_L)')
                ylabel(['Re(f_{' question5s '})'])
            
            
            elseif question6=='Imaginary'
                
                plot(p,imag(F(question5+1,:)))
                
                title(['\beta/(\gamma_L-\gamma_d)=' betais ' ' '\nu/(\gamma_L-\gamma_d)=' nuis '  ' '\alpha/(\gamma_L-\gamma_d)=' alphais '  ' '\gamma_d/\gamma_L=' gds '  ' 't\times\gamma_L=' timels])
                xlabel('s (Fourier conjugate to (kv-\omega)/\gamma_L)')
                ylabel(['Im(f_{' question5s '})'])
            
            
            end
            
            FRAMEDAT=getframe(gcf);
            aviobj = addframe(aviobj,FRAMEDAT);
        
        %Plots f_0 in real space    
        elseif question4=='Real space   '
            
            %Plots the spatially averaged part of distribution.
            plot(linspace(-Freqp/2,Freqp/2,10*npx+npx)*2*pi,real((fftshift(ifft(ifftshift([zeros(1,5*npx),F(1,:),zeros(1,5*npx)]))))*(10*npx+npx)*dp/sqrt(2*pi))+1/pi*2*pi*linspace(-Freqp/2,Freqp/2,10*npx+npx))
            % To plot it with without the equilibrium you should use:
            %plot(linspace(-Freqp/2,Freqp/2,10*npx+npx)*2*pi,real((fftshift(ifft(ifftshift([zeros(1,5*npx),F(1,:),zeros(1,5*npx)]))))*(10*npx+npx)*dp/sqrt(2*pi)))
            xlim([-10,10])
            ylim([-10/pi,10/pi])
            
            title(['\beta/(\gamma_L-\gamma_d)=' betais ' ' '\nu/(\gamma_L-\gamma_d)=' nuis '  ' '\alpha/(\gamma_L-\gamma_d)=' alphais '  ' '\gamma_d/\gamma_L=' gds '  ' 't\times\gamma_L=' timels])
            xlabel('(kv-\omega)/\gamma_L')
            ylabel('F_0 + f_0')
            
                    
        FRAMEDAT=getframe(gcf);
        aviobj = addframe(aviobj,FRAMEDAT);
            
        end
        
        
        
           else
        
        clear F2 F3 F4
     %Produces the phase space plot
        
        F_real=transpose(fftshift(ifft(ifftshift(transpose([zeros(nx,4*npx),F,zeros(nx,4*npx)])))))*(8*npx+npx)*dp/sqrt(2*pi);
        F_real=F_real(:,start:finish);
        %This adds the slope from the equilibrium distribution 
        F_real2(:,:,1)=repmat(F_real(1,:)+1/pi*O,size(x,2),1);
        for JJ=2:nx
        F_real2(:,:,JJ)=repmat(F_real(JJ,:),size(x,2),1).*(exp(i*(JJ-1)*xx'+i*pi/2*(JJ-1)));
        end
        F_real2=real(F_real2);
        %F_real2=F_real2.*M+conj(F_real2).*conj(M);
        F_real_tot=sum(F_real2,3);

             
                 
        set(gcf,'renderer','zbuffer')
        imagesc(x,O,transpose(real(F_real_tot)))
        %You can also use a contourf to make the phase plot if you really
        %want contour information, but this slows down the code a lot!
        %contourf(xx,OO,transpose(real(F_real_tot)),100,'LineStyle','none')
        set(gca,'YDir','normal')
        ylabel('(kv-\omega)/\gamma_L')
        xlabel('\xi')
        title(['t\times\gamma_L=' timels])
        try colormap(mapn)
        catch exception
        end
        %The colour bar limits should be set based on what frequency you
        %expect to be at by the end of the simulation.
        %Recall that F_0=delta_omega/pi in these units, so e.g if you have a
        %hole and you want to go to delta_omega=20, then the upper limit
        %here should be about 6.4.
        caxis([Cmin Cmax])
        colorbar
        
        FRAMEDAT=getframe(gcf);
        aviobj = addframe(aviobj,FRAMEDAT);
        
    end
    
end
        
 %Check for numerical stability  
 if p(end)*dt*abs(E2)/2>1
    break
end

end



%Closes the movie file, IMPORTANT: if the code crashes you must perform
%this close manually.  Alternatively just do CLEAR ALL
if question==1 | question==2
aviobj = close(aviobj);
end

%Saves the simulation data
save(['BOT_a=' alphais '_n=' nuis '_b=' betais '_gd=' gds '_pmax=' num2str(loopi) '_N=' num2str(nx-1) '_smax=' num2str(pbox) '_spoints=' num2str(npx) '_(' num2str(cl(1)) '_' num2str(cl(2)) '_' num2str(cl(3)) '_' num2str(cl(4)) '_' num2str(cl(5)) ').mat'],'E','F','nx','npx','pbox','loopi','gd','alphai','betai','nui','time','p','dt','dp','Freq','Freqp','nt','p0')


%Stops the clock to time the simulation
toc

%This creates the spectrogram.  This loop will continue until you click
%'cancel' on the user input screen.  This lets you create multiple
%spectrograms with different window parameters.  You can use the zoom
%features in the figure to investigate the current spectrogram before you
%click 'ok' to create a new one.

if exist('specgram') ~=0

if question==0
    options.Interpreter = 'tex';
    options.Default = 'Yes';
    question0=questdlg('Do you want to see a Fourier spectrogram','Spectrogram','Yes','No ',options);
end

while question0=='Yes'
     sig=E.*exp(-i*62.8*time);
    %2nd argument in specgram is the size of the fourier transform window (power of 2 
    %above 256), the larger this is the better resolution in frequency but 
    %the poorer in time, 4th argument should be the same as 2nd, Hann is 
    %the type of window used.  The last argument gives the window overlap
    
%Obtains spectrogram parameters from user input    
prompt={'FFT size','Window overlap'};
name='Spectrogram input';
numlines=1;
%Default values for spectrogram parameters
defaultanswer={'2048','2000'};
options.WindowStyle='normal';
answer=inputdlg(prompt,name,numlines,defaultanswer,options);
    
if isempty(answer)~=1
    figure('nextplot','add')
    specgram(real(sig),str2num(answer{1}),Freq,hann(str2num(answer{1})),str2num(answer{2}))
    title(['\beta/(\gamma_L-\gamma_d)=' betais ', ' '\nu/(\gamma_L-\gamma_d)=' nuis ', ' '\alpha/(\gamma_L-\gamma_d)=' alphais ', ' '\gamma_d/\gamma_L=' gds ', N=' nxs ', s_{max}=' pboxs ', dt\times\gamma_L=' dts ', ' npxs ' s points, dt/ds=' nts])
    xlabel('t\times \gamma_L')
    ylabel('\omega/\gamma_L/2\pi')
elseif isempty(answer)==1
    question0=0;
end

end

end



