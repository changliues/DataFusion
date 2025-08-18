BSIV_correct_model <- list(g.xz = "gn ~ x1+x2+b",
                           b.o = "b ~ x1+x2",
                           b.e = "b ~ x1+x2",
                           a.be = "a ~ x1+b",
                           a.bo = "a ~ x1+x2+b+x1:x2+x1:b+x2:b+I(x1^2)+I(x2^2)",
                           m.abo = "m ~ x1+x2+b+a+I(x1^2)+I(x2^2)+x1:x2+x1:b+x2:b+x1:a+x2:a+a:b",
                           y.abo = "y ~ x1+x2+b+a+I(x1^2)+I(x2^2)+x1:x2+x1:b+x2:b+x1:a+x2:a+a:b",
                           m.abe = "m ~ x1+x2+b+a+I(x1^2)+I(x2^2)+x1:x2+x1:b+x2:b+x1:a+x2:a+a:b")



BSIV_incorrect_model <- list(g.xb = "gn ~ x1s+x2s+b",
                             b.o = "b ~ x1s+x2s",
                             b.e = "b ~ x1s+x2s",
                             a.be = "a ~ x1s+b",
                             a.bo = "a ~ x1s+x2s+b",
                             m.abo = "m ~ x1s+x2s+b+a",
                             y.abo = "y ~ x1s+x2s+b+a",
                             m.abe = "m ~ x2s+b+a")
BSIV_ls_mod_labels <- list(
  # all correct
  mod_labels0 = c(1,1,1,1,1,1,1,1),
  # set 1 correct
  mod_labels1 = c(0,0,0,0,1,1,1,1),
  # set 2 correct
  mod_labels2 = c(1,1,1,1,1,0,0,0),
  # set 3 correct
  mod_labels3 = c(1,0,0,1,1,1,1,0),
  # set 4 correct
  mod_labels4 = c(0,1,1,0,1,0,0,1),
  # all incorrect
  mod_labels5 = c(0,0,0,0,0,0,0,0))

BSIV_model_spec <- list()

for(nlbs in 1: length(BSIV_ls_mod_labels)){
  BSIV_mod_labels <- BSIV_ls_mod_labels[[nlbs]]
  BSIV_model_spec[[nlbs]] <- list()
  for(lbs in 1:length(BSIV_mod_labels)){
    if(BSIV_mod_labels[lbs] == 1){
      BSIV_model_spec[[nlbs]][[lbs]] <- BSIV_correct_model[[lbs]]
    }
    else{
      BSIV_model_spec[[nlbs]][[lbs]] <- BSIV_incorrect_model[[lbs]]
    }
  }
  names(BSIV_model_spec[[nlbs]]) <- c("g.xb", "b.o", "b.e", "a.be", "a.bo", "m.abo", "y.abo", "m.abe")
}

BSIV_model_spec[[length(BSIV_ls_mod_labels)+1]] <- BSIV_correct_model
BSIV_model_spec[[length(BSIV_ls_mod_labels)+2]] <- BSIV_incorrect_model

names(BSIV_model_spec) <- c("all_true", "case1", "case2", "case3", "case4", "all_false", "BSIV_correct_model", " BSIV_incorrect_model")

saveRDS(BSIV_model_spec, "BSIV_model_spec.rds")
