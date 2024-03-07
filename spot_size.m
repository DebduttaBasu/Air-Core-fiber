%% Spot-Size calaulation of OAM modes 
%   Author: Debdutta Basu from IIT Madras
%   Date: 16.01.2024
clc;close all;clear all;
L=1550*(10^-9); % Wavelength
f=5*(10^-3); %lens focal length
r1=(-11:0.1:11)*(10^-6);% Checking the intensity around the focal point in transverse plane
inp_rad=0.00346/2; % aperture of lens
l=5;
w=(sqrt(l+1)*L*f)/(pi*inp_rad);
% Substituting the fixed values in the formula, parameter v dependent on f
% and r is defined

U1=exp(-(r1/w).^2).*((sqrt(2)*r1/w).^abs(l)); % Field equation at focal plane
for d1=1:length(r1) % Check for NaN values
    if isnan(U1(d1))
        U1(d1)=0;
    end
end
% Defining the x-axis, 7.38um is the newport camera cell size
%x=v*7.38*(10^-6);
x1=r1;
I1=(abs(U1.^2));
%figure,imagesc(I1);
%Normalising the intensity
y1=I1/max(max((I1)));
% Finding MFD at 1/e^2
I2=ones(length(r1));
I2=I2/((exp(1)).^2);
plot(x1*10^6,y1,'r','LineWidth',2);
hold on
%plot(x*10^6,y,'b');
%hold on
plot(x1*10^6,I2,'b','LineWidth',2);
legend('Beam Profile','1/e^2 line');%'Method-2(Lommel Function)',
% Increase line width and font size of the current axis (gca)
xlabel('x (micrometers)');
ylabel('Intensity');
ax = gca;
ax.LineWidth = 1; % Set the line width
ax.FontSize = 20; % Set the font size
grid on;
