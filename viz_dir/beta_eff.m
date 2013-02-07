function k = beta_eff(Te)
k  = (Te<=1200) .* (1.95e-7*(300./Te).^.7) + (Te>1200) .* (1.73e-7*( 300./Te).^.61); % cm6/s, Lillis, 2009
end