function I_t2 = JadduShimemura1_I_t2(t,t1,y1_t1)
%JADDUSHIMEMURA1_I_T2
%    I_T2 = JADDUSHIMEMURA1_I_T2(T,T1,Y1_T1)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    22-Oct-2019 20:05:45

I_t2 = (t.*1.6e+1+(t-1.0./2.0).^2.*8.0-1.7e+1./2.0).^2./2.0e+2+((t-1.0./2.0).^2.*8.0-1.0./2.0).^2+(t.*(3.0./2.0)-t1.*(3.0./2.0)+y1_t1-t.^2.*4.0+t.^3.*(8.0./3.0)+t1.^2.*4.0-t1.^3.*(8.0./3.0)).^2;