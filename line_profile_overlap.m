%% code for overlap integral between free space LG modes with OAM beam inside the air-core fiber
%   Author: Debdutta Basu from IIT Madras
%   Date: 16.01.2024
clc;clear all;close all;
%% Conversion to cylindrical cordinates
xpixel=600; % Defining the number of pixels along x
ypixel=600; % Defining the number of pixels along y
cell=0.05e-6; % Cell size of each pixel
x=-(xpixel/2):1:(xpixel/2); x=x*cell;
y=-(ypixel/2):1:(ypixel/2); y=y*cell;
[X,Y]=meshgrid(x,y);
rho = sqrt(Y.^2+X.^2);    % [phi,rho]=cart2pol(X,Y) is an equivalent function
phi =atan2(-Y,X);        % -Y is used to care of the reverse number mapping
phi(phi<0)=phi(phi<0)+2*pi; % This command is for mapping the phase profile from 0 to 2pi
%E_field_even_y=0.028494298563127.*sin(phi);
load('All_mode.mat');
%% Loading HE even and odd modes
HE_even_x= mode19_Ex;            %load("Mode_4_Ex.csv");
HE_even_y= mode19_Ey;            %load("Mode_4_Ey.csv");
HE_odd_x=  mode20_Ex;            %load("Mode_5_Ex.csv");
HE_odd_y=  mode20_Ey;            %load("Mode_5_Ey.csv");
% figure;imshow(imag(HE_even_x));colormap('turbo');c=colorbar;c.FontSize = 18;
% title('HE_{6,1}^{Odd}','FontSize',18,'FontWeight', 'bold');
% figure;imshow(abs(HE_odd_x));colormap('turbo');c=colorbar;c.FontSize = 18;
% title('HE_{6,1}^{Even}','FontSize',18,'FontWeight', 'bold');
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');

oam_x=HE_even_x+1i.*HE_odd_x;
%oam_y=1*(-1i*HE_even_y+HE_odd_y);
oam_y=-1.*(HE_even_y+1i.*HE_odd_y);
norm_x= sqrt(sum(sum(abs(oam_x).^2)));
norm_oam_x = oam_x/norm_x;
norm_y= sqrt(sum(sum(abs(oam_y).^2)));
norm_oam_y = oam_y/norm_y;
oam=(1/sqrt(2))*(oam_x+oam_y);
norm_1= sqrt(sum(sum(abs(oam).^2)));
norm_oam = oam/norm_1;
%oam = oam/max(max(oam));
% figure;imagesc(abs(norm_oam)/max(max(abs(norm_oam))));colormap('turbo'),colorbar;
% xlabel('Pixel');
% ylabel('Pixel');
% figure;imshow(abs(norm_oam)/max(max(abs(norm_oam))));colormap('turbo');c=colorbar;c.FontSize = 18;
% title('Intensity of OAM_{+5}','FontSize',18,'FontWeight', 'bold');
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% figure;imagesc(abs(oam_x));olormap('turbo'),colorbar;
% xlabel('Pixel');
% ylabel('Pixel');
% figure;imagesc(angle(oam_x));
% figure;imshow(angle(oam_x));colormap('jet');c=colorbar;c.FontSize = 18;
% title('Phase of OAM_{+5}','FontSize',18,'FontWeight', 'bold');
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% figure;imagesc(abs(oam_x));colormap('turbo'),colorbar;
% figure;imagesc(angle(oam_y));
% figure;imagesc(angle(oam_y));
oam_intensity= abs(oam.^2);%sqrt(abs(oam_x).^2+abs(oam_y).^2);
% figure;imagesc(abs(oam));
% figure;imagesc(oam_intensity)
% figure;
% plot(x*1e6,oam_intensity(:,300),'r','LineWidth',2);

%% Defining the modes : LG modes
L=1550*(10^-9); % Wavelength
f=5*(10^-3); %lens focal length
inp_rad=0.00346/2; % aperture of lens
l=5;
w=(sqrt(l+1)*L*f)/(pi*inp_rad);
  % input beam which is a plane wave; we are converting it into beam with waist 0.8mm 
p=0; %radial order =0,Normalisation constant M is substituted l=1,p=0 and put in equation
cal_field_re=(sqrt(2/pi))*(1/w).*((sqrt(2)*rho/w).^abs(l)).*exp(-(rho/w).^2).*exp(1i*(l)*phi);
norm_cal_field_re= cal_field_re/sqrt(sum(sum(abs(cal_field_re).^2)));
I1=(abs(norm_cal_field_re.^2));
y1=I1/max(max((I1)));
% figure;imagesc(abs(norm_cal_field_re)/max(max(abs(norm_cal_field_re))));colormap('turbo'),colorbar;
% xlabel('Pixel');
% ylabel('Pixel');
% figure;
% plot(x*1e6,y1(:,300),'r','LineWidth',2);
%% Overlap Integral
E1 =norm_oam;
E2 =norm_cal_field_re;
% Define the integrand function as the dot product of E1 and E2
integrand = abs((sum(sum(E1.*conj(E2))))^2);
norm_factor = 1;%sum(sum((E1)^2))*sum(sum((E2)^2));
OI= (integrand/norm_factor)
f_values = 5*(10^-3);%linspace(4e-3, 8e-3, 10); % Variation of focal length from 4 to 8

for idx = 1:length(f_values)
    f = f_values(idx); % Current focal length
    w = (sqrt(l + 1) * L * f) / (pi * inp_rad);
    cal_field_re = (sqrt(2 / pi)) * (1 / w) .* ((sqrt(2) * rho / w).^abs(l)) .* exp(-(rho / w).^2) .* exp(1i * l * phi);
    norm_cal_field_re = cal_field_re / sqrt(sum(sum(abs(cal_field_re).^2)));    
    I1 = (abs(norm_cal_field_re.^2));
    y1 = I1 / max(max((I1)));
    % figure;
    % plot(x * 1e6, y1(:, 300), 'r', 'LineWidth', 2);
    
    %% Overlap Integral
    E1 = norm_oam;
    E2 = norm_cal_field_re;
    integrand = abs((sum(sum(E1 .* conj(E2))))^2);
    norm_factor = sum(sum((E1).^2)) * sum(sum((E2).^2));
    OI = (integrand / norm_factor);
    overlap_integral_values(idx) = OI;
end









% Perform the triple integral using the integral3 function
%overlap_integral = integral2(integrand,min(min(x)),max(max(x)),min(min(y)),max(max(y)));

%disp(['Overlap Integral: ', num2str(overlap_integral)]);
line=ones(601,601);
line=line/((exp(1)).^2);
figure;
plot(x,abs((oam_x(300,:))),'r','LineWidth',2);
hold on
plot(x,abs((oam_y(300,:))),'b','LineWidth',2)
hold on 
plot(x,abs((line(300,:))),'k','LineWidth',2);
% Increase line width and font size of the current axis (gca)
xlabel('x (micrometers)');
ylabel('Intensity');
ax = gca;
ax.LineWidth = 1; % Set the line width
ax.FontSize = 14; % Set the font size
% 
% %% Calculating the orthonormalisation coefficient
% Oam_1 = (oam_x.*conj(oam_x))+(oam_y.*conj(oam_y));
% Oam_2 = (fftshift(fft2((Oam_1))));
% Oam_3 = Oam_2*(cell)^2; % cell is the conversion factor
% figure; imagesc(abs(Oam_3));title('Coefficients')
% oam_fft = abs(Oam_3(251,251));
% Efield_OAM_Beam=oam/(sqrt(oam_fft));
% Efield_OAM_Beam_real=oam_x/(sqrt(oam_fft));
% Efield_OAM_Beam_im=oam_y/(sqrt(oam_fft));
% Efield_OAM_Beam_total_1=Efield_OAM_Beam_real+Efield_OAM_Beam_im;
% Efield_OAM_Beam_total=sqrt(abs(Efield_OAM_Beam_real).^2+abs(Efield_OAM_Beam_im).^2);
% figure; imagesc(abs(Efield_OAM_Beam_total))
% % Efield_OAM_Beam_intensity=sqrt(((abs(oam_x).^2+abs(oam_y).^2))./oam_fft);
% % figure;imagesc(angle(Efield_OAM_Beam));title('Normalised OAM beam phase');
% % figure;imagesc(Efield_OAM_Beam_intensity);title('Normalised OAM beam Intensity');
% 
% %((((sqrt(2/pi))*(1/w).*((sqrt(2)*rho/w).^abs(l)).*exp(-(rho/w).^2).*exp(1i*(l)*phi))+(exp(1i*pi/2)*((sqrt(2/pi))*(1/w).*((sqrt(2)*rho/w).^abs(l)).*exp(-(rho/w).^2).*exp(1i*(l)*phi))))/(sqrt(LG_Coefficient_M2_fft))).*conj(((((sqrt(2/pi))*(1/w).*((sqrt(2)*rho/w).^abs(l)).*exp(-(rho/w).^2).*exp(1i*(l)*phi))+(exp(1i*pi/2)*((sqrt(2/pi))*(1/w).*((sqrt(2)*rho/w).^abs(l)).*exp(-(rho/w).^2).*exp(1i*(l)*phi))))/(sqrt(LG_Coefficient_M2_fft))))
% %% Method 2 : using FFT to calculate the normalization coefficient
% %cal_field_LG=cal_field.*conj(cal_field);
% EF_LG = (fftshift(fft2((cal_field_LG))));
% EF_LG=EF_LG*(cell)^2; % cell is the conversion factor
% %figure; imagesc(abs(EF_LG));title('Coefficients')
% LG_Coefficient_M2_fft=abs(EF_LG(251,251));
% Efield_LG_Beam=cal_field/(sqrt(LG_Coefficient_M2_fft));
% Efield_LG_Beam_real=cal_field_re/(sqrt(LG_Coefficient_M2_fft));
% Efield_LG_Beam_im=cal_field_im/(sqrt(LG_Coefficient_M2_fft));
% Efield_LG_Beam_total_1=Efield_LG_Beam_real+Efield_LG_Beam_im;
% Efield_LG_Beam_total=sqrt(abs(Efield_LG_Beam_real).^2+abs(Efield_LG_Beam_im).^2);
% figure;imagesc((Efield_LG_Beam_total));title('Normalised LGbeam intensity');
% figure;
% plot(abs((Efield_LG_Beam_total(250,:))),'r')
% hold on
% plot(abs((Efield_OAM_Beam_total(250,:))),'b')
% return;
% Efield_LG_Beam_in=sqrt(((abs(cal_field_re).^2+abs(cal_field_im).^2))./LG_Coefficient_M2_fft);
% %figure;imagesc(angle(Efield_LG_Beam));title('Normalised LGbeam phase');
% %figure;imagesc(abs(Efield_LG_Beam_in));title('Normalised LGbeam intensity');
% %% Performing overlap
% overlap_LG_oam_1=((cal_field_re.*conj(oam_x))./(((sqrt(LG_Coefficient_M2_fft)))*((sqrt(oam_fft)))))+((cal_field_im.*conj(oam_y))./(((sqrt(LG_Coefficient_M2_fft)))*((sqrt(oam_fft)))));
% overlap_LG_oam_2=(fftshift(fft2((overlap_LG_oam_1))));
% overlap_LG_oam_3=overlap_LG_oam_2*(cell)^2; % cell is the conversion factor
% %figure; imagesc(abs(overlap_LG_oam_3));title('Coefficients')
% overlap_LG_oam=abs(overlap_LG_oam_3(251,251));
% a=8*10^-6;
