
# test that evenness decreases when cv_abund increases

library(mobsim)
library(dplyr)
library(purrr)
library(ggplot2)


n_total<-1e3
iter = 1e3
ev_test<- map_dfr(seq(30, 60, 5), function(rich){
  map_dfr(10^seq(-1, 2, 0.2), function(cva){
    map_dfr(1:iter, function(iter){
      ab = as.numeric(sim_sad(s_pool = rich, n_sim =n_total
                              , sad_type = "lnorm", sad_coef = list(cv_abund = cva)))
      ev = (MeanRarity::rarity(ab, l = -1)-1)/(MeanRarity::rarity(ab, l = 1) -1)
      return(data.frame( rich, cva, ev))
    })
  })
})

pdf(file = "cv_abund_and_evenness.pdf")
ev_test %>% ggplot(aes(cva, ev, color  = rich))+
  geom_jitter(alpha = 0.3)+
  geom_smooth(alpha = 0.5, aes(group = rich), se = F) +
  theme_classic() +
  labs(x = "cv_abund", color = "s_pool", y = "evenness" ) +
  scale_x_log10()
dev.off()

# yes but super noisy
