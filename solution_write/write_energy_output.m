function write_energy_output(PRO,CON,Wint_n,Wext_n,WKE,Wvdamp_n,time)
    
string='a';

if (~PRO.rest && CON.incrm==1)
    string='w';
    system('rm energy.dat');
    system('rm VDenergy.dat');
end

string3=sprintf('energy.dat');

fid5= fopen(string3,string);

%--------------------------------------------------------------------------
% Print energy information.
%--------------------------------------------------------------------------
fprintf(fid5,'%5.5e %5.5e  %5.5e %5.5e %5.5e\n',time, WKE, Wint_n, Wext_n,Wvdamp_n);

fclose(fid5);

% string4=sprintf('VDenergy.dat');
% fid4= fopen(string4,string);
% fprintf(fid5,'%5.5e %5.5e\n',time, Wvdamp_n);
% fclose(fid4);
end