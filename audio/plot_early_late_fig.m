function plot_early_late_fig()


        %%
        figure
        subplot(3,2,1)
        x = delta_pos_rise*bin_size_up;
        y = 'pos rise time';
        histogram(x(2:end))
        if ~isnan(x(1))
            xline(x(1),'r');
            sig = inv_prcnt(x);
            title({[y ': late - early'],['p=' num2str(sig)]})
        else
            title([y ': late - early'])
        end
        xlabel('\Delta(m)')
        
        subplot(3,2,2)
        x = delta_neg_rise*bin_size_up;
        y = 'neg rise time';
        histogram(x(2:end))
        if ~isnan(x(1))
            xline(x(1),'r');
            sig = inv_prcnt(x);
            title({[y ': late - early'],['p=' num2str(sig)]})
        else
            title([y ': late - early'])
        end
        xlabel('\Delta(m)')
        
        subplot(3,2,3)
        x = delta_pos_sig*bin_size;
        y = 'pos 1st sig bin';
        histogram(x(2:end))
        if ~isnan(x(1))
            xline(x(1),'r');
            sig = inv_prcnt(x);
            title({[y ': late - early'],['p=' num2str(sig)]})
        else
            title([y ': late - early'])
        end
        xlabel('\Delta(m)')
        
        subplot(3,2,4)
        x = delta_neg_sig*bin_size;
        y = 'neg 1st sig bin';
        histogram(x(2:end))
        if ~isnan(x(1))
            xline(x(1),'r');
            sig = inv_prcnt(x);
            title({[y ': late - early'],['p=' num2str(sig)]})
        else
            title([y ': late - early'])
        end
        xlabel('\Delta(m)')
        
        subplot(3,2,5)
        x = delta_peak*bin_size;
        y = 'peak ego tuning';
        histogram(x(2:end))
        if ~isnan(x(1))
            xline(x(1),'r');
            sig = inv_prcnt(x);
            title({[y ': late - early'],['p=' num2str(sig)]})
        else
            title([y ': late - early'])
        end
        xlabel('\Delta(m)')
        
        %%
        %         figure
        %         subplot(1,2,1)
        %         hold on
        %         plot(dis_X_bins_vector_of_centers,ego_early.ego_shuffle(2:end,:),'color',[0.9 0.9 0.9],'LineWidth',3);
        %         plot(dis_X_bins_vector_of_centers,ego_early.ego_shuffle(1,:),'color',[1 0 1],'LineWidth',3);
        %         %             plot(ego_early.x_bins_up,ego_early.mid_shuf_up,'color',[0 0 0],'LineWidth',1.5);
        %         if ~isempty(ego_early.sig_pos_bins)
        %             plot(dis_X_bins_vector_of_centers(ego_early.sig_pos_bins),ego_early.ego_shuffle(1,ego_early.sig_pos_bins),'g*')
        %         end
        %         if ~isempty(ego_early.sig_neg_bins)
        %             plot(dis_X_bins_vector_of_centers(ego_early.sig_neg_bins),ego_early.ego_shuffle(1,ego_early.sig_neg_bins),'r*')
        %         end
        %         if ~isempty(ego_early.width_pos)
        %             plot(ego_early.x_bins_up(ego_early.width_pos),ego_early.ego_true_up(ego_early.width_pos),'b')
        %         end
        %         if ~isempty(ego_early.width_neg)
        %             plot(ego_early.x_bins_up(ego_early.width_neg),ego_early.ego_true_up(ego_early.width_neg),'b')
        %         end
        %         title('early')
        %         subplot(1,2,2)
        %         hold on
        %         plot(dis_X_bins_vector_of_centers,ego_late.ego_shuffle(2:end,:),'color',[0.9 0.9 0.9],'LineWidth',3);
        %         plot(dis_X_bins_vector_of_centers,ego_late.ego_shuffle(1,:),'color',[1 0 1],'LineWidth',3);
        %         %             plot(ego_late.x_bins_up,ego_late.mid_shuf_up,'color',[0 0 0],'LineWidth',1.5);
        %         if ~isempty(ego_late.sig_pos_bins)
        %             plot(dis_X_bins_vector_of_centers(ego_late.sig_pos_bins),ego_late.ego_shuffle(1,ego_late.sig_pos_bins),'g*')
        %         end
        %         if ~isempty(ego_late.sig_neg_bins)
        %             plot(dis_X_bins_vector_of_centers(ego_late.sig_neg_bins),ego_late.ego_shuffle(1,ego_late.sig_neg_bins),'r*')
        %         end
        %         title('late')
        %         if ~isempty(ego_late.width_pos)
        %             plot(ego_late.x_bins_up(ego_late.width_pos),ego_late.ego_true_up(ego_late.width_pos),'b')
        %         end
        %         if ~isempty(ego_late.width_neg)
        %             plot(ego_late.x_bins_up(ego_late.width_neg),ego_late.ego_true_up(ego_late.width_neg),'b')
        %         end
        
        
        %%
        %         figure
        %         p=panel();
        %         sgtitle([cell_name ', dir_' num2str(ii_dir)],'interpreter','none')
        %         p.pack({0.2,[]},2);
        %         %     p.select('all')
        %         %     p.identify();
        %         p(1,1).select()
        %         title('early')
        %         plot(firing_rate.dis_X_bins_vector_of_centers{1,1},dis_x_pos_fr_for_2D_early,'m','LineWidth',2)
        %         p(2,1).select()
        %         fn_plot_2D_field(allo_ego_map_early,firing_rate(ii_dir).dis_X_bins_vector{1},firing_rate(ii_dir).dis_X_bins_vector_of_centers{1},firing_rate(ii_dir).allo_X_bins_vector{1},firing_rate(ii_dir).allo_X_bins_vector_of_centers{1},[]);
        %         axis xy
        %         p(1,2).select()
        %         title('late')
        %         plot(firing_rate.dis_X_bins_vector_of_centers{1,1},dis_x_pos_fr_for_2D_late,'m','LineWidth',2)
        %         p(2,2).select()
        %         fn_plot_2D_field(allo_ego_map_late,firing_rate(ii_dir).dis_X_bins_vector{1},firing_rate(ii_dir).dis_X_bins_vector_of_centers{1},firing_rate(ii_dir).allo_X_bins_vector{1},firing_rate(ii_dir).allo_X_bins_vector_of_centers{1},[]);
        %         axis xy
        %         p.margin = [15 15 15 20];