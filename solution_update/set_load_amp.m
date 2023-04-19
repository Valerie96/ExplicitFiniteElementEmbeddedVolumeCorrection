function CON = set_load_amp(CON,tn)
len = size(CON.ramp_table,1);

for i=1:len-1
    if tn>=CON.ramp_table(i,1) && tn<=CON.ramp_table(i+1,1)
        t0 = CON.ramp_table(i,1); t1=CON.ramp_table(i+1,1);
        A0 = CON.ramp_table(i,2); A1=CON.ramp_table(i+1,2);
    end
end

An = (tn-t0)/(t1-t0)*(A1-A0)+A0;
CON.amp=An;