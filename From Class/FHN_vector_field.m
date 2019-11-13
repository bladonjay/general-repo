function FHN_vector_field(v,w, I,a,b,tau)

  v0 = (min(v):0.01:max(v));
  w_for_v_nullcline = v0 - 1/3 * v0.^3 + I;
  w_for_w_nullcline = (v0-a)./b;
  
  plot(v0,w_for_v_nullcline,'--')
  hold on
  plot(v0,w_for_w_nullcline,'--')

  [v0,w0] = meshgrid(v, w);

  dv = v0 - 1/3.*v0.^3 - w0 + I;
  dw = (v0 - a - b.*w0)/tau;
  
  quiver(v0,w0,dv,dw)
  hold off
  xlabel('v')
  ylabel('w')
  xlim([min(v) max(v)])
  ylim([min(w) max(w)])
  
end