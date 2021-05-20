joint_model_fit_1_new <- joint_model_fit_1
load(file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_wrkng/newfit.RData'))
load(file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_wrkng/oldfit.RData'))

all.equal(joint_model_fit_1, joint_model_fit_1_new)

all.equal(joint_model_fit_1$mcmc$bs_gammas$`1`, joint_model_fit_1_new$mcmc$bs_gammas$`1`)
