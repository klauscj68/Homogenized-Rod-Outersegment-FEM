% Conversion factors for longitudinal diffusion and pointwise diffusion
% coefficients

flag_data='mouse';
switch flag_data
    case 'salamander'

        % Salamander data
        R=5.5;
        sigma=15/14;
        epsilon_0=0.014;
        H=22.4;
        A_inc=0.8;
        theta=1/2;
        l=epsilon_0/theta;
        A_gap=2*pi*R*sigma*epsilon_0;

    case 'mouse'

        % Mouse data
        R=0.7;
        sigma=15/14;
        epsilon_0=0.014;
        H=23.6;
        A_inc=0.2593*0.3111/2;
        theta=1/2;
        l=epsilon_0/theta;
        A_gap=2*pi*R*sigma*epsilon_0;

end 



% Holcman and Korenbrot formula
gamma_h=(A_inc+A_gap)*l/(pi*R^2*l/2+(A_inc+A_gap)*l/2);



% Our formula (also of Lamb)
f_A=(A_inc+A_gap)/(pi*R^2+A_gap+A_inc);
f_V=((1-theta)*pi*R^2*H+A_inc*H+A_gap*H)/(pi*R^2*H);
gamma_o=f_A/f_V;