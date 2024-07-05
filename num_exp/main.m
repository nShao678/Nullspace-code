% clear all

%% Table 5.3
% N = 7;
% nd = 10;
% result = cell(1,N);
% for ii = 1:N
%     tic
%     eps0 = 4^(-3-ii);
%     result{ii} = main_GL_SV(eps0,nd);
%     toc
% end
% save('data_GL_SV')
% table = zeros(4,7);
% for ii = 1:7
%     table(1,ii) = result{ii}.eps0;
%     table(2,ii) = result{ii}.lambdaN;
%     table(3,ii) = result{ii}.normAV;
%     table(4,ii) = result{ii}.normVAV;
% end
% latex(sym(table));

%% Fig 5.1
% N = 7;
% nd = 10;
% eps0 = 4^(-10);
% result = main_GL_SVH(eps0,nd);
% xx = 1:result.zero(end);
% yy = zeros(size(xx));
% yy(end) = result.nz;
% for ii = 1:result.nz-1
%     idx = result.zero(ii):result.zero(ii+1)-1;
%     yy(idx) = ii;
% end
% figure
% hold on
% plot(xx,yy,'-','LineWidth',2)
% ylabel('Number of eigenvalues of $T_\ell$ below $3\epsilon$','Interpreter','latex')
% xlabel('Dimension $\ell$ of Krylov subspace','Interpreter','latex')
% set(gcf, 'Color', 'w');
% hold off
% export_fig('conHis.eps')
% export_fig('conHis.pdf')
% save('data_GL_SVH')

%% Table 5.4
% nameSet = {'GL7d12','Franz7','Franz9','Franz10','GL7d13','Franz8'};
% check = [128,256,1024,1024,2048,3072];
% 
% nBmax = [5,6,8,8,9,9];
% iiMax = length(nameSet);
% result = cell(1,iiMax);
% for ii = 1:iiMax
%     name = ['SSmatrix/',nameSet{ii}]
%     result{ii} = main_CH(name,check(ii),nBmax(ii));
% end
% save('data_CH')
% for ii = 1:6
%     len = length(result{ii}.matvec);
%     idx = len-4:len;
%     latex(sym([result{ii}.matvec(idx),0;result{ii}.time(idx),result{ii}.nulltime]))
% end


%% table 5.5

% result = main_GL_LB(3);
% save('data_GL_LB')
% 
% latex(sym([result.matvec;result.time;result.err]))


% result = main_CH_LB;
% save('data_CH_LB')
% latex(sym([result.matvec;result.time;result.err]))


%% table 5.6

% result = main_CH_14;
% save('data_CH_14')
% latex(sym([result.matvec;result.time]))




