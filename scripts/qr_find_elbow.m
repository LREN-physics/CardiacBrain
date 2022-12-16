function ElbowPoint=qr_find_elbow(x,y,PlotFig)
% This function find the elbow point of a L-Curve using the “Kneedle” algorithm
%
% ElbowPoint=qr_find_elbow(x,y,PlotFig)
% INPUTS:
%     - x is the abscissa of the L-curve
%     - y is the L-curve
% Optional :
%     - PlotFig is a 0/1 flag to show a figure or not
% OUTPUTS:
%     - ElbowPoint is the value of the elbow along x
%
% The code has been addapted from Satopää, Albrecht, Irwin, and Raghavan (2011)
% https://towardsdatascience.com/detecting-knee-elbow-points-in-a-graph-d13fc517a63c
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

if nargin<3
    PlotFig=0;
end

%% Putting everything in a 1D vector
x=x(:);
y=y(:);

Line = interp1([x(1),x(end)],[y(1),y(end)],x);

%% Normalization
nx=(x-min(x))./max(x-min(x));
ny=(y-min(y))./max(y-min(y));
nLine=(Line-min(Line))./max(Line-min(Line));

%% Calculating Segment1 :
% The L2 distance between curve and droite
nDroiteHighSampling=linspace(nLine(1),nLine(end),1000);
nxHighSampling=linspace(0,1,1000);

Segment1=zeros(size(nLine));
for cx=1:length(nLine)
    temp=zeros(length(nDroiteHighSampling),1);
    for cy=1:length(nDroiteHighSampling) % L2 normal calculation can be improved
        temp(cy)=norm(nxHighSampling(cy)+1i*nDroiteHighSampling(cy)-(nx(cx)+1i*ny(cx)));
    end
    Segment1(cx)=min(temp(:));
end

%% Calculating Segment2
% The difference between curve and droite
Segment2=nLine-ny;

%% Elbow point
% Finding the position of the maximum difference between difference and L2
% norm, which is the elbow point
[~,nxmax]=max(Segment2-Segment1);
ElbowPoint=x(nxmax);

%% Figures
if PlotFig==1
%     figure
%     hold on
%     plot(nx,Segment1)
%     plot(nx,Segment2)
%     plot([nx(nxmax),nx(nxmax)],[0 1])

    figure
    hold on
    plot(nx,ny,'.-')
    plot(nx,nLine,'--')
    plot(nx,Segment2-Segment1,'.-')
    plot([nx(nxmax),nx(nxmax)],[0 1])
%     daspect([1 1 1])
    legend('Normalized curve','Given line','Difference curve',['Elbow point = ',num2str(nx(nxmax))])
    set(gcf,'color','w');
    set(gca,'fontsize',14)
end

end