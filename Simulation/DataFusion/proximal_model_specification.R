## labels indicating whether the following models are correctly or incorrectly specified (in order)
# 1. g.x.fit
# 2. a.ex.fit
# 3. eta.fit


proximal_correct_model <- list(g.x = "g ~ x1+x2+x3",
                               a.ex = "ae ~ x1e+x3e",
                               eta = "w ~ x10e+x20e+x30e+I(x10e^2)+I(x20e^2)+x10e:x20e+x10e:x30e+x30e:x20e+I(x10e^3)+I(x20e^3)+I(x10e^2):x20e+I(x10e^2):x30e+I(x20e^2):x10e+I(x20e^2):x30e+x10e:x20e:x30e")



proximal_incorrect_model <- list(g.x = "g ~ x1s+x3s",
                                 a.ex = "ae ~ x1se+x2se",
                                 eta = "w ~ x1s0e+x2s0e+x3s0e")

proximal_ls_mod_labels <- list(
  # all correct
  mod_labels0 = c(1,1,1),
  # set 1 correct
  mod_labels1 = c(0,0,1),
  # set 2 correct
  mod_labels2 = c(1,1,0),
  # set 3 correct
  mod_labels3 = c(1,1,0),
  # all incorrect
  mod_labels4 = c(0,0,0))

proximal_model_spec <- list()

for(nlbs in 1: length(proximal_ls_mod_labels)){
  proximal_mod_labels <- proximal_ls_mod_labels[[nlbs]]
  proximal_model_spec[[nlbs]] <- list()
  for(lbs in 1:length(proximal_mod_labels)){
    if(proximal_mod_labels[lbs] == 1){
      proximal_model_spec[[nlbs]][[lbs]] <- proximal_correct_model[[lbs]]
    }
    else{
      proximal_model_spec[[nlbs]][[lbs]] <- proximal_incorrect_model[[lbs]]
    }
  }
  names(proximal_model_spec[[nlbs]]) <- c("g.x", "a.ex", "eta")
}

proximal_model_spec[[length(proximal_ls_mod_labels)+1]] <- proximal_correct_model
proximal_model_spec[[length(proximal_ls_mod_labels)+2]] <- proximal_incorrect_model

names(proximal_model_spec) <- c("all_true", "case1", "case2", "case3", "all_false", "proximal_correct_model", "proximal_incorrect_model")

saveRDS(proximal_model_spec, "proximal_model_spec.rds")
