function CloseFigures(fig_nums)

for ii=1:length(fig_nums)
   close(figure(fig_nums(ii)))    
end

end

